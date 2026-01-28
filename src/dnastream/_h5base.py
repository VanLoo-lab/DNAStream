from abc import ABC, abstractmethod
import h5py
import warnings
import pandas as pd
import numpy as np
from .schema import Schema
from .utils import as_str
from typing import Callable, Any

from .constants import (
    Hook,
)


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

        self._schema = None
        self._hooks: list[Hook] = []

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
    def schema(self) -> Schema | None:
        return self._schema

    @property
    def path(self) -> str:
        """str: Absolute HDF5 path to the dataset."""
        return f"{self._parent.name}/{self._name}"

    def register_hook(self, hook: Hook) -> None:
        self._hooks.append(hook)

    def _emit(self, scope: str, event: str, fn: str, **payload: Any) -> None:
        for hook in self._hooks:
            hook(scope=scope, event=event, dataset=self.path, fn=fn, **payload)

    def exists(self) -> bool:
        """Check whether the dataset exists in the parent group.

        Returns
        -------
        bool
            True if the dataset exists on disk, otherwise False.
        """
        return self._name in self._parent

    def open(
        self, schema: Schema | None = None, *, strict: bool = True
    ) -> h5py.Dataset:
        """Open the underlying HDF5 dataset.

        Validates schema identity against metadata stored on disk.

        Parameters
        ----------
        schema : dnastream.Schema or None, optional
            Schema used to validate schema identity. If None, uses the cached schema
            from the last successful `create()` or `open()`.
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

        if schema is None:
            if self._schema is None:
                raise ValueError(
                    f"No cached schema for '{self.name}'. Pass a Schema to open(), or call create(schema=...)."
                )
            schema = self._schema

        if not isinstance(schema, Schema):
            raise ValueError("schema must be of type 'dnastream.Schema'")
        if not self.exists():
            raise RuntimeError(
                f"Dataset does not exist at {self.path}. Cannot be opened!"
            )

        self._validate_schema(schema, strict)
        # Cache schema for downstream calls (e.g., add/update can call open() with schema=None).
        self._schema = schema

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

    def _validate_schema(self, expected_schema: Schema, strict: bool) -> None:
        """Validate dataset schema identity against an expected schema.

        Parameters
        ----------
        expected_schema : Schema
            Expected schema containing ``schema_version`` and
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

        ds_schema_version = as_str(ds_schema_version)
        ds_schema_hash = as_str(ds_schema_hash)

        expected_schema_version = as_str(expected_schema.version)
        expected_schema_hash = as_str(expected_schema.hash())

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
        schema: Schema,
        *,
        if_exists: str = "raise",
        validate: bool = True,
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
        schema : Schema
            dnastream.Schema object
        if_exists : {"raise", "open"}, optional
            the behavior if the dataset already exists
        validation : bool, optional
            if the dataset should be validated
        shape : tuple[int, ...], default (0,)
            Initial dataset shape.
        maxshape : tuple[int | None, ...], default (None,)
            Maximum shape for resizable datasets.
        **kwargs
            Forwarded to `h5py.Group.create_dataset`, e.g. `shuffle=True`, etc.

        Raises
        ------
        RuntimeError
            If the dataset already exists and ``if_exists=raise``.
        ValueError
            If ``if_exists`` is invalid.

        Returns
        -------
        h5py.Dataset
            The created dataset.
        """

        # can't create a DNAStream dataset without a schema
        if not isinstance(schema, Schema):
            raise ValueError("schema must be of type 'dnastream.Schema'")

        if if_exists in ["raise", "open"]:

            if self.exists():
                if if_exists == "raise":
                    raise RuntimeError(
                        f"Refusing to create: {self.path} already exists."
                    )
                else:
                    if validate:
                        self.validate()
                    warnings.warn(
                        f"Registry aleady exists at {self.path}, returning existing handle.",
                        stacklevel=2,
                    )
                    self._schema = schema
                    return self.open(schema)

        else:
            raise ValueError(
                f"if_exists {if_exists} not valid, expected one of ['raise', 'open']"
            )

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

        ds = self._parent.create_dataset(
            self._name,
            shape=shape,
            maxshape=maxshape,
            dtype=schema.dtype,
            compression=compression,
            compression_opts=compression_opts,
            chunks=chunks,
            **kwargs,
        )

        # Persist schema identity on disk
        ds.attrs["schema_version"] = str(schema.version)
        ds.attrs["schema_hash"] = schema.hash()
        ds.attrs["schema_pairs_json"] = schema.json_pairs()

        ds.attrs["name"] = str(self._name)

        self._schema = schema

        if validate:
            self.validate()
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

    def to_dataframe(self, arr: np.ndarray | None = None, **_: Any) -> pd.DataFrame:
        if arr is None:
            arr = self._ds()[:]
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
            idx = np.asarray(selector, dtype=object)
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
        df = self.to_dataframe(arr=arr)

        # restore original requested order (optional but usually expected)
        inv = np.empty_like(order)
        inv[order] = np.arange(order.size)
        return df.iloc[inv].reset_index(drop=True)
