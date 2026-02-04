import h5py
import os
import getpass
import socket
from datetime import datetime, timezone
import uuid
from .registry import Registry
from ._builtin_schemas import REGISTRY_SCHEMAS, PROVENANCE_SCHEMAS
from .provenance import Provenance
from .io import IO
import logging
from typing import Literal

from .constants import Event, SCOPE, EVENTS

from .utils import _qualname, package_version


class DNAStream:
    def __init__(self, path, verbose=False):
        """
        The core class for a DNAStream object.

        :param path: str
            the path to the HDF5 file
        :param mode: str (one of 'r', 'r+', 'x')
        """

        self.path = path

        self._mode = None

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
        schema = attrs.get("package_version", None)
        creator = attrs.get("created_by", None)
        creator_host = attrs.get("created_on_host", None)
        patient = attrs.get("patient_id", None)
        created_at = attrs.get("created_at", None)

        mystr = f"{fmt} for Patient {patient} with Package Version: {schema}\n"
        mystr += f"File path: '{self.path}'\n"
        mystr += f"Created by: {creator} on host {creator_host} at {created_at}\n"
        mystr += f"Mode: '{self.mode}'"
        mystr += "\nDataset groups:\n"
        for table, ref in self.handle.items():
            mystr += f"\t'{table}': {len(ref)} tables\n"

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

    def set_patient_id(self, patient_id):
        """
        Update the patient_id in the HDF5 file.

        patient_id : str
            The new patient_id to write into the DNAStream file.
        """
        old_patient_id = self.handle.attrs["patient_id"]
        self.handle.attrs["patient_id"] = str(patient_id)
        self._emit(
            EVENTS.MODIFY,
            _qualname(self.set_patient_id),
            old_patient_id=old_patient_id,
            patient_id=patient_id,
        )

    def registry(self, name: str) -> Registry:
        """Returns a reference to a Registry object by name.

        name : str
        The name of the requested registry.
        """
        if name not in self._registries:
            raise KeyError(f"Registry with key '{name}' does not exist.")
        return self._registries[name]

    def provenance(self, name: str) -> Provenance:
        """Returns a reference to a Provenance object by name.

        name : str
        The name of the requested provenance object.
        """
        if name not in self._provenance:
            raise KeyError(f"Proveance with key '{name}' does not exist.")
        return self._provenance[name]

    def _bind_components(self, strict=True):
        self._load_provenance(strict=True)
        self._load_registries(strict=True)
        self._register_hooks()

    @classmethod
    def create(cls, path: str, *, patient_id: str = "", verbose: bool = False):
        """Create and initialize a new DNAStream HDF5 file.

        path : str
            The path to where the DNAStream file should be created.
        patient_id : str, optional
            The patient id attribute for the DNAStream file
        verbose : bool, optional
            Flag to provide more verbose output.


        Raises
        ------
        FileExistsError : if a file exists at the provided path
        """
        if os.path.isfile(path):
            raise FileExistsError(
                f"'{path} already exists, DNAStream file not created, use `open()` instead."
            )
        ds = cls(path, verbose)
        ds._initialize_new_file(patient_id=patient_id)

        ds._bind_components(strict=True)

        ds._emit(EVENTS.CREATE, _qualname(ds.create), file=path)

        if ds.verbose:
            ds.logger.info(f"Created DNAStream file {path}")
        return ds

    @classmethod
    def open(
        cls,
        path: str,
        *,
        mode: Literal["r", "r+"] = "r",
        verbose: bool = False,
    ):
        """Open the stream to DNAStream file.

        path : str
            The path to where the DNAStream file should be created.
        mode : str, optional
            The mode used to open the stream.  'r' is read-only and 'r+' gives read/write access.
        verbose : bool, optional
            Flag to provide more verbose output.

        Raises
        ------
        FileNotFoundError : if file does not exist at path

        ValueError : if an invalid mode is passed
        """

        if mode not in ("r", "r+"):
            raise ValueError(
                f"Mode '{mode}' invalid for opening stream; use 'r' or 'r+'."
            )
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f"File '{path}' not found. Use `create()` to initialize file before opening."
            )

        ds = cls(path, verbose=verbose)

        try:
            ds._handle = h5py.File(name=path, mode=mode)
            ds._mode = mode

            ds._validate_header()
            ds._bind_components(strict=True)
            return ds

        except Exception:
            # avoid leaking a live handle on failure
            ds.close()
            raise

    def is_connected(self) -> bool:
        """Returns a boolean indicating whether or not the stream is connected.."""
        return (
            self._handle is not None
            and getattr(self._handle, "id", None) is not None
            and self._handle.id.valid
        )

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
            raise FileNotFoundError(
                f"DNAStream file '{path}' does not exist. Use 'create()' to initialize it."
            )

    @property
    def handle(self):
        """Return the live DNAStream HDF5 handle; raise if not connected."""
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

    def __enter__(self):
        """Enter a context where the underlying HDF5 handle is open.

        Notes
        -----
        - For existing files, prefer modes 'r' or 'r+'.
        - For new files, call `create()` (mode 'x'/'w-') before using the context manager.

        This method calls `connect()` and then verifies that a live handle exists.
        """

    def __enter__(self):
        if not self.is_connected():
            raise RuntimeError(
                "Use ds.open(mode=...) or ds.create(...) before entering context."
            )
        if self.verbose:
            self.logger.info(f"Connection to DNAStream file '{self.path}' open.")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _load_provenance(self, *, strict: bool = True) -> None:
        """Bind Provenance wrappers for provenance log

        For new files (mode 'x'), this will create any missing provenance datasets.
        For existing files (modes 'r'/'r+'), this requires the datasets to exist.
        """
        # if self.mode in ("x", "w-"):
        #     prov_grp = self.handle.require_group("provenance")

        if "provenance" not in self.handle:
            raise ValueError("DNAStream file is missing required group '/provenance'.")
        prov_grp = self.handle["provenance"]

        for key, schema in PROVENANCE_SCHEMAS.items():
            prov = Provenance(prov_grp, key)
            if key in prov_grp:
                # Existing dataset: open + validate schema identity
                prov.open(schema, strict=strict)

            self._provenance[key] = prov

    def _load_registries(self, *, strict: bool = True) -> None:
        """Bind Registry wrappers for all known registries.

        For new files (mode 'x'), this will create any missing registry datasets.
        For existing files (modes 'r'/'r+'), this requires the datasets to exist.
        """

        reg_grp = self.handle["registry"]

        for key, schema in REGISTRY_SCHEMAS.items():
            reg = Registry(reg_grp, key)

            if key in reg_grp:
                # Existing dataset: open + validate schema identity
                reg.open(schema, strict=strict)

            self._registries[key] = reg

    def _initialize_new_file(self, patient_id="") -> None:
        """Initialize the DNAStream layout on an already-open, empty HDF5 handle."""
        # Safety: only initialize empty files/handles

        self._handle = h5py.File(self.path, mode="x")

        if len(self.handle.keys()) != 0 or len(self.handle.attrs) != 0:
            raise RuntimeError("Refusing to initialize: file is not empty.")

        # Required file-level attributes
        self.handle.attrs["dnastream_format"] = "DNAStream"
        # self.handle.attrs["schema_version"] = PACKAGE_VERSION
        self.handle.attrs["package_version"] = package_version("dnastream")
        self.handle.attrs["created_by"] = getpass.getuser()
        self.handle.attrs["created_on_host"] = socket.gethostname()
        self.handle.attrs["created_at"] = (
            datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        )

        # Recommended
        self.handle.attrs["file_uuid"] = str(uuid.uuid4())
        # self.handle.attrs["software_version"] = (
        #     PACKAGE_VERSION  # or package __version__
        # )

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

        prov_grp = self.handle["provenance"]
        for key, schema in PROVENANCE_SCHEMAS.items():
            prov = Provenance(prov_grp, key)
            prov.create(schema)

        self._load_provenance()

        reg_grp = self.handle["registry"]
        for key, schema in REGISTRY_SCHEMAS.items():
            reg = Registry(reg_grp, key)
            reg.create(schema)

            self._record_event(
                SCOPE.DNASTREAM,
                EVENTS.CREATE,
                dataset=reg.path,
                fn=_qualname(reg.create),
                schema_version=schema.version,
                schema_hash=schema.hash(),
            )
        self._load_registries()

    @property
    def patient_id(self):
        return self.handle.attrs["patient_id"]

    @property
    def mode(self):
        return self._mode

    def _validate_header(self):
        """Validate that this is DNAStream file before modification or reading."""
        fmt = self.handle.attrs.get("dnastream_format", None)

        if fmt != "DNAStream":
            raise ValueError(
                "File is not specified as a DNAStream object: missing/invalid attribute 'dnastream_format'."
            )
        if "package_version" not in self.handle.attrs:
            raise ValueError(
                "DNAStream file missing required attribute 'schema_version'."
            )

    def _handle_provider(self):
        return self._handle if self.is_connected() else None

    def _record_event(self, scope: str, event: Event, dataset: str, fn: str, **payload):
        self.provenance("log").add(scope, event, dataset, fn=fn, **payload)

    def _emit(self, event, fn, **payload):
        self._record_event(SCOPE.DNASTREAM, event, ".", fn, **payload)

    def _register_hooks(self):
        for _, reg in self._registries.items():
            reg.register_hook(self._record_event)
