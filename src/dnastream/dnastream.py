import h5py
import os
import getpass
import socket
from datetime import datetime, timezone
import uuid


SCHEMA_VERSION = "0.1.0"
PACKAGE_VERSION = "0.1.0"


class DNAStream:
    def __init__(self, path, mode, patient=""):
        """
        The core class for a DNAStream object.

        :param path: str
            the path to the HDF5 file
        :param mode: str (one of 'r', 'r+', 'x')
        """
        self._validate_mode(mode)
        self.path = path
        self.mode = mode
        self._handle = None
        self._connected = False
        self.patient = patient

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

    @staticmethod
    def _validate_mode(mode):
        """
        Validate that the specified mode is a valid DNAStream mode.

        Allowed modes:
            - 'r'  : read-only
            - 'r+' : read/write existing
            - 'x'  : create new (fail if exists)

        Mode 'w' is intentionally disallowed to prevent silent data loss.
        """
        if mode not in ("r", "r+", "x"):
            raise ValueError(
                "Requested operating mode does not exist.\n"
                "DNAStream does not support mode 'w' because it can silently overwrite data.\n"
                "Use:\n"
                "  'x'  to create a new file (fails if file exists)\n"
                "  'r+' to modify an existing file\n"
                "  'r'  to read only"
            )

    @staticmethod
    def _validate_path(path):
        """
        Validate that the h5 file path exists

        :param path: str
            path to the h5 file
        """
        if not os.path.isfile(path):
            raise FileNotFoundError(f"HDF5 file '{path}' does not exist.")

    def connect(self):
        if (
            self._handle is not None
            and getattr(self._handle, "id", None) is not None
            and self._handle.id.valid
        ):
            return

        if self.mode in ("r", "r+"):
            self._validate_path(self.path)

        self._handle = h5py.File(self.path, self.mode)

        if self.mode == "x":
            self._initialize_new_file()
        else:
            self._validate_header()

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

    def __enter__(self):
        self.connect()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _initialize_new_file(self):
        """Initialize a new DNAStream file (mode 'x')."""

        # Safety: only initialize empty files/handles
        if len(self.handle.keys()) != 0 or len(self.handle.attrs) != 0:
            raise RuntimeError("Refusing to initialize: file is not empty.")

        # Required file-level attributes
        self.handle.attrs["dnastream_format"] = "DNAStream"
        self.handle.attrs["schema_version"] = SCHEMA_VERSION
        self.handle.attrs["created_by"] = getpass.getuser()
        self.handle.attrs["created_on_host"] = socket.gethostname()
        self.handle.attrs["created_at"] = (
            datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        )

        # Recommended
        self.handle.attrs["file_uuid"] = str(uuid.uuid4())
        self.handle.attrs["dnastream_software_version"] = (
            PACKAGE_VERSION  # or package __version__
        )

        # Domain
        self.handle.attrs["patient_id"] = self.patient

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

        # Reserve provenance subgroups
        self.handle["provenance"].require_group("runs")
        self.handle["provenance"].require_group("changes")

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
