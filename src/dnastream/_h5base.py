from abc import ABC, abstractmethod
import h5py
import json
import warnings
import pandas as pd
import numpy as np


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
        if "/" in name or name == "":
            raise ValueError("Dataset name must be a non-empty relative name (no '/').")
        self._parent = parent
        self._name = name
        self._validated = False

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
            If True, schema mismatches and validation errors (with validate=True) raise ``ValueError``. If False, mismatches
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

        if not self._validated:
            # Validate subclass-specific invariants
            try:
                self.validate()

            except Exception as e:
                if strict:
                    raise
                warnings.warn(f"Validation failed for {self.path}: {e}", stacklevel=2)

        self._validated = True
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
        self,
        schema: dict,
        *,
        shape=(0,),
        maxshape=(None,),
        compression="gzip",
        compression_opts=4,
        chunks=True,
        **kwargs,
    ) -> h5py.Dataset:
        """
        Create the underlying HDF5 dataset under this object's parent group.

        Parameters
        ----------
        schema : dict
            Compiled schema dict. Must include key "dtype".
        shape : tuple[int, ...], default (0,)
            Initial dataset shape.
        maxshape : tuple[int | None, ...], default (None,)
            Maximum shape for resizable datasets.
        **kwargs
            Forwarded to `h5py.Group.create_dataset`, e.g. `shuffle=True`, etc.

        Returns
        -------
        h5py.Dataset
            The created dataset.
        """
        if self.exists():
            raise RuntimeError(f"Refusing to create: {self.path} already exists.")

        # Prevent callers from overriding core invariants
        forbidden = {
            "name",
            "dtype",
            "shape",
            "maxshape",
            "compression",
            "compression_opts",
            "chunks",
        }
        overlap = forbidden.intersection(kwargs)
        if overlap:
            raise TypeError(
                f"Do not pass {sorted(overlap)} to create(); "
                f"they are controlled by H5Dataset."
            )

        if "dtype" not in schema:
            raise ValueError("schema must contain key 'dtype'")

        ds = self._parent.create_dataset(
            self._name,
            shape=shape,
            maxshape=maxshape,
            dtype=schema["dtype"],
            compression=compression,
            compression_opts=compression_opts,
            chunks=chunks,
            **kwargs,
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

    @abstractmethod
    def validate(self, *args, **kwargs):
        """Validates the contract of the subclass.

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

    @staticmethod
    def _to_dataframe(arr: np.ndarray) -> pd.DataFrame:
        # Decode fixed-width byte string fields (dtype kind 'S') before constructing the DataFrame.
        if isinstance(arr, np.ndarray) and arr.dtype.names is not None:
            arr2 = arr.copy()
            for name in arr2.dtype.names:
                dt = arr2.dtype[name]
                if dt.kind == "S":
                    arr2[name] = np.char.decode(arr2[name], "utf-8")
            return pd.DataFrame(arr2)

        # Fallback for non-structured arrays
        return pd.DataFrame(arr)

    def get(
        self, selector=None, *, by: str | None = None, mode: str = "active_only"
    ) -> pd.DataFrame:
        mode = self._validate_mode(mode)
        ds = self._ds()

        if selector is None:
            # still avoid ds[:] if you want: implement to_dataframe(mode) similarly with on-disk mask
            return self.to_dataframe(mode=mode)

        if callable(selector):
            # this path inherently needs a DataFrame of *something* to call the predicate on
            # in the future consider to disallow callables for lazy mode
            df = self.to_dataframe(mode=mode)
            return df[selector(df)]

        if by is None:
            raise ValueError(
                "Must pass by='label', by='id', or by='idx' when selector is not callable."
            )

        if by == "label":
            idx = self.resolve_idx_from_labels(selector, mode=mode)
        elif by == "id":
            idx = self.resolve_idx_from_ids(selector, mode=mode)
        elif by == "idx":
            idx = np.asarray(selector)
        else:
            raise ValueError("by must be one of {'label','id','idx'}")

        # normalize scalar -> 1D
        idx = np.asarray(idx, dtype=object)
        if idx.ndim == 0:
            idx = np.array([idx], dtype=object)

        # drop missing and convert to int
        idx = idx[~pd.isna(idx)].astype(int)

        order = np.argsort(idx)
        idx_sorted = idx[order]

        arr = ds[idx_sorted]
        df = self._to_dataframe(arr)

        # restore original requested order (optional but usually expected)
        inv = np.empty_like(order)
        inv[order] = np.arange(order.size)
        return df.iloc[inv].reset_index(drop=True)
