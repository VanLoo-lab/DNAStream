from abc import ABC, abstractmethod
import h5py
import json
import warnings


class H5Dataset(ABC):
    """Abstract base class for lazily accessed HDF5 datasets.

    ``H5Dataset`` represents a *named dataset within an existing HDF5 group*.
    It provides shared logic for:

    - Checking dataset existence
    - Creating datasets with a declared schema
    - Opening datasets with optional schema validation
    - Persisting schema identity metadata on disk

    It stores a reference to the parent :class:`h5py.Group` and the
    dataset name, and resolves the dataset lazily when accessed.

    Subclasses must implement :meth:`add` to define how new records are appended
    to the dataset.

    Parameters
    ----------
    parent : h5py.Group
        Parent HDF5 group in which the dataset resides.
    name : str
        Name of the dataset within ``parent``.

    Attributes
    ----------
    group : h5py.Group
        The parent group containing this dataset.
    group_name : str
        Absolute HDF5 path of the parent group.
    name : str
        Dataset name relative to the parent group.
    path : str
        Absolute HDF5 path to the dataset.

    Notes
    -----
    - Construction does not require that the dataset already exists on disk.
    - Schema identity is tracked using HDF5 attributes (``schema_version`` and
      ``schema_hash``).
    - Subclasses implement the concrete append semantics in :meth:`add`.

    See Also
    --------
    h5py.Group
    h5py.Dataset
    """

    def __init__(self, parent: h5py.Group, name: str) -> None:
        """Construct an H5Dataset bound to a parent group and dataset name.

        Parameters
        ----------
        parent : h5py.Group
            Parent HDF5 group that will contain the dataset.
        name : str
            Dataset name relative to ``parent``.

        Notes
        -----
        Construction does not create the dataset on disk. Use :meth:`create` to
        create it, and :meth:`open` to access it.
        """
        self._parent = parent
        self._name = name

    def __str__(self) -> str:
        """Return a human-readable summary of the dataset."""
        if not self.exists():
            return f"{self.path} (missing)"
        ds = self.open()
        return f"{self.path} (Shape: {ds.shape})"

    @property
    def group(self) -> h5py.Group:
        """h5py.Group: Parent group that contains this dataset."""
        return self._parent

    @property
    def group_name(self) -> str:
        """str: Absolute HDF5 path of the parent group."""
        return self._parent.name

    @property
    def name(self) -> str:
        """str: Dataset name relative to the parent group."""
        return self._name

    @property
    def path(self) -> str:
        """str: Absolute HDF5 path to the dataset."""
        return f"{self._parent.name}/{self._name}"

    def exists(self) -> bool:
        """Check whether the dataset exists in the parent group.

        Returns
        -------
        bool
            True if the dataset exists on disk, otherwise False.
        """
        return self._name in self._parent

    def open(
        self, expected_schema: dict | None = None, *, strict: bool = True
    ) -> h5py.Dataset:
        """Open the underlying HDF5 dataset.

        Optionally validates schema identity against metadata stored on disk.

        Parameters
        ----------
        expected_schema : dict or None, optional
            Compiled schema dictionary (e.g., from ``dnastream.schemas.compile_schema``)
            containing at least ``schema_version`` and ``schema_hash``. If provided,
            these values are compared against attributes stored on the dataset.
        strict : bool, optional
            If True, schema mismatches raise ``ValueError``. If False, mismatches
            emit a warning and the dataset is still returned.

        Returns
        -------
        h5py.Dataset
            The opened dataset handle.

        Raises
        ------
        RuntimeError
            If the dataset does not exist.
        ValueError
            If schema validation fails and ``strict=True``.
        """
        if not self.exists():
            raise RuntimeError(
                f"Dataset does not exist at {self.path}. Cannot be opened!"
            )
        if expected_schema is not None:
            self._validate_schema(expected_schema, strict)
        return self._parent[self._name]

    def _ds(self) -> h5py.Dataset:
        """Return the underlying dataset without performing validation.

        Returns
        -------
        h5py.Dataset
            Dataset handle.

        Notes
        -----
        This is an internal helper. Call :meth:`open` for user-facing access and
        optional schema validation.
        """
        return self._parent[self._name]

    def _validate_schema(self, expected_schema: dict, strict: bool) -> None:
        """Validate dataset schema identity against an expected schema.

        Parameters
        ----------
        expected_schema : dict
            Expected schema dictionary containing ``schema_version`` and
            ``schema_hash`` keys.
        strict : bool
            If True, mismatches raise ``ValueError``. If False, mismatches emit a
            warning.

        Raises
        ------
        ValueError
            If schema attributes are missing or mismatched and ``strict=True``.
        """
        ds = self._ds()
        ds_schema_version = ds.attrs.get("schema_version", None)
        ds_schema_hash = ds.attrs.get("schema_hash", None)

        ds_schema_version = (
            None if ds_schema_version is None else str(ds_schema_version)
        )
        ds_schema_hash = None if ds_schema_hash is None else str(ds_schema_hash)

        expected_schema_version = expected_schema.get("schema_version", None)
        expected_schema_version = (
            None if expected_schema_version is None else str(expected_schema_version)
        )
        expected_schema_hash = expected_schema.get("schema_hash", None)
        expected_schema_hash = (
            None if expected_schema_hash is None else str(expected_schema_hash)
        )

        problems: list[str] = []

        missing = False
        if not ds_schema_version or not expected_schema_version:
            problems.append("missing schema attr schema_version on dataset/schema")
            missing = True
        if not ds_schema_hash or not expected_schema_hash:
            problems.append("missing schema attr schema_hash on dataset/schema")
            missing = True

        if not missing:
            if expected_schema_hash != ds_schema_hash:
                problems.append(
                    f"schema_hash mismatch: expected {expected_schema_hash[:12]}…, found {ds_schema_hash[:12]}…"
                )
            if expected_schema_version != ds_schema_version:
                problems.append(
                    f"schema_version mismatch: expected {expected_schema_version}, found {ds_schema_version}"
                )

        if problems:
            msg = f"Schema check failed for {self.path}: " + "; ".join(problems)
            if strict:
                raise ValueError(msg)
            warnings.warn(msg)

    def create(
        self, schema: dict, *, shape=(0,), maxshape=(None,), chunks=True
    ) -> h5py.Dataset:
        """Create the dataset on disk and persist schema identity attributes.

        Parameters
        ----------
        schema : dict
            Compiled schema dictionary. Must include a ``dtype`` key (a NumPy dtype
            suitable for an HDF5 dataset). If present, ``schema_version``,
            ``schema_hash``, and ``schema_pairs`` are stored as HDF5 attributes.
        shape : tuple, optional
            Initial dataset shape (default is ``(0,)``).
        maxshape : tuple, optional
            Maximum dataset shape (default is ``(None,)``).
        chunks : bool or tuple, optional
            Chunking strategy passed to :meth:`h5py.Group.create_dataset`.

        Returns
        -------
        h5py.Dataset
            The newly created dataset handle.

        Raises
        ------
        RuntimeError
            If the dataset already exists.
        KeyError
            If ``schema['dtype']`` is missing.
        """
        if self.exists():
            raise RuntimeError(f"Refusing to create: {self.path} already exists.")

        ds = self._parent.create_dataset(
            self._name,
            shape=shape,
            maxshape=maxshape,
            dtype=schema["dtype"],
            chunks=chunks,
        )

        # Persist schema identity on disk
        if "schema_version" in schema:
            ds.attrs["schema_version"] = str(schema["schema_version"])
        if "schema_hash" in schema:
            ds.attrs["schema_hash"] = str(schema["schema_hash"])
        if "schema_pairs" in schema:
            ds.attrs["schema_pairs_json"] = json.dumps(
                list(schema["schema_pairs"]),
                separators=(",", ":"),
            )

        ds.attrs["name"] = str(self._name)

        return ds

    @abstractmethod
    def add(self, *args, **kwargs):
        """Append records to the dataset.

        Subclasses must implement this method to define how new data are
        validated, transformed, and appended to the underlying HDF5 dataset.

        Returns
        -------
        Any
            Implementation-defined. Common choices are ``None`` or an
            :class:`h5py.Dataset` handle after appending.

        Notes
        -----
        - Implementations should assume the dataset already exists (i.e.,
          :meth:`create` has been called at least once).
        - This base class does not prescribe the input format; concrete dataset
          types (e.g., registries, measurements, provenance tables) should
          document their accepted inputs.
        """
        ...
