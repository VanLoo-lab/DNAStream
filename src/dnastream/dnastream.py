import h5py
import os
import getpass
import warnings
import socket
from datetime import datetime, timezone
import uuid
from .registry import Registry
from ._builtin_schemas import REGISTRY_SCHEMAS, PROVENANCE_SCHEMAS
from .provenance import Provenance
from .io import IO
import logging
from typing import Literal, Callable, Any

from .constants import PACKAGE_VERSION, Event, SCOPE, EVENTS

from .utils import _qualname


class DNAStream:
    def __init__(self, path, mode, verbose=False):
        """
        The core class for a DNAStream object.

        :param path: str
            the path to the HDF5 file
        :param mode: str (one of 'r', 'r+', 'x')
        """

        self.path = path

        mode = self._validate_mode(mode)

        self.mode = mode

        self._handle = None

        self._registries = {}
        self._provenance = {}

        self.io = IO(self)
        self.verbose = verbose

        logging.basicConfig(
            format="{asctime} - {levelname} - {message}",
            style="{",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.INFO,
        )
        self.logger = logging.getLogger()

    def __str__(self):
        if self._handle is None:
            return f"DNAStream(path={self.path!r}, mode={self.mode!r}, connected=False)"

        attrs = self.handle.attrs

        fmt = attrs.get("dnastream_format", None)
        schema = attrs.get("schema_version", None)
        creator = attrs.get("created_by", None)
        creator_host = attrs.get("created_on_host", None)
        patient = attrs.get("patient_id", None)
        created_at = attrs.get("created_at", None)

        mystr = f"{fmt} for Patient {patient} with Schema Version: {schema}\n"
        mystr += f"File path: {self.path}\n"
        mystr += f"Created by: {creator} on host {creator_host} at {created_at}\n"
        mystr += f"Mode: {self.mode}"

        return mystr

    def __getattr__(self, name: str):
        # Allow access like ds.sample / ds.variant

        def raise_if_not_connected():
            if not self.is_connected():
                raise RuntimeError(
                    "DNAStream is not connected. Call connect() before accessing registries."
                )

        if name in REGISTRY_SCHEMAS:
            raise_if_not_connected()
            return self.registry(name)

        if name in PROVENANCE_SCHEMAS:
            raise_if_not_connected()
            return self.provenance(name)

        raise AttributeError(name)

    def set_mode(self, mode: str):
        mode = self._validate_mode(mode)
        self.mode = mode

    def create(self, patient_id="") -> None:
        """Create and initialize a new DNAStream HDF5 file.

        Requires mode 'x' or 'w-' (fail if the file already exists).
        Leaves the instance connected.
        """
        if self.mode not in ("x", "w-"):
            raise ValueError(
                f"Mode {self.mode} is not valid for creating a new DNAStream HDF5 file. Use mode 'x' or 'w-'.",
            )

        if os.path.isfile(self.path):
            warnings.warn(
                f"A valid DNAStream file already exists at {self.path}, use `connect() instead.`",
                UserWarning,
                stacklevel=2,
            )
            return
            # self.connect()

        if self.is_connected():
            return

        # Create the file and initialize the DNAStream layout
        self._handle = h5py.File(name=self.path, mode=self.mode)
        self._initialize_new_file(patient_id=patient_id)

        # Create required datasets for registries/provenance
        self._load_provenance(strict=True)
        self._load_registries(strict=True)
        self.register_hooks()

        # Now that provenance exists, log file creation
        self._emit(EVENTS.CREATE, _qualname(self.create), file=self.path)

        if self.verbose:
            self.logger.info(f"Created DNAStream file {self.path}")

    def open(self) -> None:
        """Open an existing DNAStream file (alias for connect)."""
        self.connect()

    def connect(self) -> None:
        if self.is_connected():
            return

        if self.mode in ("x", "w-"):
            warnings.warn(
                "DNAStream.connect() is for opening existing files. "
                "Mode is 'x'/'w-'; calling create() instead.",
                UserWarning,
                stacklevel=2,
            )
            self.create()
            return

        if self.mode not in ("r", "r+"):
            warnings.warn(
                f"Mode {self.mode!r} is not valid for connect(); expected 'r' or 'r+'.",
                UserWarning,
                stacklevel=2,
            )
            return

        # normal open-existing path
        self._validate_path(self.path)
        self._handle = h5py.File(name=self.path, mode=self.mode)
        self._validate_header()
        self._load_provenance(strict=True)
        self._load_registries(strict=True)
        self.register_hooks()

    @staticmethod
    def _validate_mode(mode) -> str:
        """
        Validate that the specified mode is a valid DNAStream mode.

        Allowed modes:
            - 'r'  : read-only
            - 'r+' : read/write existing
            - 'x'  : create new (fail if exists)
            - 'w-' : create new (fail if exists)

        Mode 'w' and 'a' are intentionally disallowed to prevent silent data loss.
        """
        if mode not in ("r", "r+", "x", "w-"):
            raise ValueError(
                "Requested operating mode does not exist.\n"
                "DNAStream does not support mode 'w' because it can silently overwrite data.\n"
                "Use:\n"
                "  'x/w-'  to create a new file (fails if file exists)\n"
                "  'r+' to modify an existing file\n"
                "  'r'  to read only"
            )
        return mode

    @staticmethod
    def _validate_path(path):
        """
        Validate that the h5 file path exists

        :param path: str
            path to the h5 file
        """
        if not os.path.isfile(path):
            raise FileNotFoundError(f"HDF5 file '{path}' does not exist.")

    @property
    def handle(self):
        """Return the live HDF5 handle; raise if not connected."""
        if self._handle is None:
            raise RuntimeError(
                "DNAStream is not connected. Call connect() or use it as a context manager."
            )
        return self._handle

    def close(self):
        """Close the underlying HDF5 handle."""
        if self._handle is not None:
            try:
                self._handle.close()
            finally:
                self._handle = None

        if self.verbose:
            self.logger.info(f"Connection to DNAStream file '{self.path}' closed.")

    def is_connected(self) -> bool:
        return (
            self._handle is not None
            and getattr(self._handle, "id", None) is not None
            and self._handle.id.valid
        )

    def __enter__(self):
        """Enter a context where the underlying HDF5 handle is open.

        Notes
        -----
        - For existing files, prefer modes 'r' or 'r+'.
        - For new files, call `create()` (mode 'x'/'w-') before using the context manager.

        This method calls `connect()` and then verifies that a live handle exists.
        """
        self.connect()
        if not self.is_connected():
            raise RuntimeError(
                "DNAStream context manager could not open a file handle. "
                "If you intended to create a new file, call create() first (mode 'x'/'w-'). "
                "If you intended to open an existing file, use mode 'r' or 'r+'."
            )
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _load_provenance(self, *, strict: bool = True) -> None:
        """Bind Provenance wrappers for provenance log

        For new files (mode 'x'), this will create any missing provenance datasets.
        For existing files (modes 'r'/'r+'), this requires the datasets to exist.
        """
        if self.mode in ("x", "w-"):
            prov_grp = self.handle.require_group("provenance")
        else:
            if "provenance" not in self.handle:
                raise ValueError(
                    "DNAStream file is missing required group '/provenance'."
                )
            prov_grp = self.handle["provenance"]

        for key, schema in PROVENANCE_SCHEMAS.items():
            prov = Provenance(prov_grp, key)
            created = False
            if key in prov_grp:
                # Existing dataset: open + validate schema identity
                prov.open(schema, strict=strict)
            else:
                # Missing dataset: only allowed when creating a new file
                if self.mode in ("x", "w-"):
                    prov.create(schema)
                    created = True
                else:
                    raise ValueError(
                        f"DNAStream file is missing required provenance '/proveance/{key}'."
                    )

            self._provenance[key] = prov
            # Only log creation when we actually created the dataset (new files)
            if self.mode in ("x", "w-") and created:
                self.log_event(
                    SCOPE.DNASTREAM,
                    EVENTS.CREATE,
                    dataset=prov.path,
                    fn=_qualname(prov.create),
                    schema_version=schema.version,
                    schema_hash=schema.hash(),
                )

    def _load_registries(self, *, strict: bool = True) -> None:
        """Bind Registry wrappers for all known registries.

        For new files (mode 'x'), this will create any missing registry datasets.
        For existing files (modes 'r'/'r+'), this requires the datasets to exist.
        """
        if self.mode in ("x", "w-"):
            reg_grp = self.handle.require_group("registry")
        else:
            if "registry" not in self.handle:
                raise ValueError(
                    "DNAStream file is missing required group '/registry'."
                )
            reg_grp = self.handle["registry"]

        for key, schema in REGISTRY_SCHEMAS.items():
            reg = Registry(reg_grp, key)

            if key in reg_grp:
                # Existing dataset: open + validate schema identity
                reg.open(schema, strict=strict)
            else:
                # Missing dataset: only allowed when creating a new file
                if self.mode in ("x", "w-"):
                    reg.create(schema)
                else:
                    raise ValueError(
                        f"DNAStream file is missing required registry '/registry/{key}'."
                    )

            self._registries[key] = reg

    def _initialize_new_file(self, patient_id="") -> None:
        """Initialize the DNAStream layout on an already-open, empty HDF5 handle."""
        # Safety: only initialize empty files/handles
        if len(self.handle.keys()) != 0 or len(self.handle.attrs) != 0:
            raise RuntimeError("Refusing to initialize: file is not empty.")

        # Required file-level attributes
        self.handle.attrs["dnastream_format"] = "DNAStream"
        self.handle.attrs["schema_version"] = PACKAGE_VERSION
        # self.handle.attrs["package_version"] = PACKAGE_VERSION
        self.handle.attrs["created_by"] = getpass.getuser()
        self.handle.attrs["created_on_host"] = socket.gethostname()
        self.handle.attrs["created_at"] = (
            datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        )

        # Recommended
        self.handle.attrs["file_uuid"] = str(uuid.uuid4())
        self.handle.attrs["software_version"] = (
            PACKAGE_VERSION  # or package __version__
        )

        # Domain
        self.handle.attrs["patient_id"] = patient_id

        # Reserved top-level groups
        for grp in (
            "registry",
            "measurements",
            "links",
            "results",
            "canonical",
            "provenance",
        ):
            self.handle.require_group(grp)

        # # Reserve provenance subgroups
        # self.handle["provenance"].require_group("runs")
        # self.handle["provenance"].require_group("changes")

    def set_patient_id(self, patient_id):
        old_patient_id = self.handle.attrs["patient_id"]
        self.handle.attrs["patient_id"] = str(patient_id)
        self._emit(
            EVENTS.MODIFY,
            _qualname(self.set_patient_id),
            old_patient_id=old_patient_id,
            patient_id=patient_id,
        )

    @property
    def patient_id(self):
        return self.handle.attrs["patient_id"]

    def _validate_header(self):
        """Validate that this is DNAStream file before modification or reading."""
        fmt = self.handle.attrs.get("dnastream_format", None)
        if fmt != "DNAStream":
            raise ValueError(
                "File is not specified as a DNAStream object: missing/invalid attribute 'dnastream_format'."
            )
        if "schema_version" not in self.handle.attrs:
            raise ValueError(
                "DNAStream file missing required attribute 'schema_version'."
            )

    def _handle_provider(self):
        return self._handle if self.is_connected() else None

    def registry(self, name: str) -> Registry:
        if name not in self._registries:
            raise KeyError(f"Registry with key '{name}' does not exist.")
        return self._registries[name]

    def provenance(self, name: str) -> Registry:
        if name not in self._provenance:
            raise KeyError(f"Proveance with key '{name}' does not exist.")
        return self._provenance[name]

    def log_event(self, scope: str, event: Event, dataset: str, fn: str, **payload):
        self.provenance("log").add(scope, event, dataset, fn=fn, **payload)

    def _emit(self, event, fn, **payload):
        self.log_event(SCOPE.DNASTREAM, event, ".", fn, **payload)

    def register_hooks(self):
        for _, reg in self._registries.items():
            reg.register_hook(self.log_event)
