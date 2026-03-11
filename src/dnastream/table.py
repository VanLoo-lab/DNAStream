"""dnastream.table

This module implements an HDF5-backed append-only 1-D table

The primary entry point is :class:`~dnastream.table.Table`.
"""

from __future__ import annotations
from typing import Callable
import uuid
import getpass
from datetime import datetime, timezone
import numpy as np
import h5py

from ._h5base import H5Dataset
from .utils import norm_key, as_str, decode_arr

import pandas as pd
from .schema import Schema
from .constants import (


    EVENTS,
)


class Table(H5Dataset):
    """HDF5-backed append-only table.

    A `Table` stores rows in an HDF5 compound dataset. Each row has a unique
    identifier (`id`). 


    Parameters
    ----------
    parent : h5py.Group
        Parent HDF5 group.
    name : str
        Dataset name.



    """

    def __init__(self, parent: h5py.Group, name: str):
        """Create a new `Table` wrapper.

        Parameters
        ----------
        parent : h5py.Group
            Parent group.
        name : str
            Dataset name.

        Notes
        -----
        This constructor does not create the dataset. It only binds to an HDF5 parent.
        """
 
        super().__init__(parent, name)

   

    def __len__(self):
        """Number of rows in the table."""
        return self._ds().shape[0]

    def __repr__(self) -> str:
        n = "?"  # avoid disk access if dataset isn't available
        try:
            n = len(self)
        except Exception:
            pass
        return f"{self.__class__}(name={self.name!r}, path={self.path!r}, " f"n={n})"

    def __iter__(self):
        for row in self._ds():
            yield decode_arr(row)

    def __contains__(self, item) -> bool:
        if item is None:
            return False

        item = norm_key(as_str(item))

        ds = self.open()
        names = ds.dtype.names or ()
        if "id" not in names:
            raise RuntimeError(f"'id' field not present in table. Table '{self.name}' is invalid")

        if ds.shape[0] == 0:
            return False

        ids = ds["id"][:]                    
        ids_norm = {norm_key(x) for x in ids}  # normalize bytes/str consistently
        return item in ids_norm


    def __getitem__(self, item) -> dict[str, object]:
        """Return a table row (as a dict) by id.

        Parameters
        ----------
        item
            Row id (UUID, str, or bytes). Must be a scalar.

        Returns
        -------
        dict[str, object]
            Mapping from field name to value. Byte fields are decoded as UTF-8 and
            NumPy scalars are converted to Python scalars.

        Raises
        ------
        KeyError
            If the id is not present.
        TypeError
            If `item` is not a scalar id-like value.
        RuntimeError
            If the table does not contain an 'id' field.
        """
        if item is None:
            raise KeyError("None is not a valid id.")

        # Unwrap numpy 0-d scalars
        if isinstance(item, np.ndarray) and item.shape == ():
            item = item.item()

        # Accept common scalar id-like types
        if not isinstance(item, (str, bytes, np.bytes_, np.str_, uuid.UUID)):
            raise TypeError(
                f"Id lookup expects a scalar str/bytes/UUID; got {type(item).__name__}."
            )

        rid = norm_key(as_str(item))

        ds = self.open()
        names = ds.dtype.names or ()
        if "id" not in names:
            raise RuntimeError(
                f"'id' field not present in table. Table '{self.name}' is invalid"
            )

        if ds.shape[0] == 0:
            raise KeyError(f"id {rid!r} not found in table '{self.name}'.")

        # Cache-free: scan ids on disk (materialize the id column once)
        ids = ds["id"][:]
        idx = None
        for i, v in enumerate(ids.tolist()):
            if norm_key(v) == rid:
                idx = i
                break

        if idx is None:
            raise KeyError(f"id {rid!r} not found in table '{self.name}'.")

        return decode_arr(ds[idx])

    @property
    def columns(self) -> tuple[str, ...]:
        ds = self.open()
        return tuple(ds.dtype.names or ())

    @property
    def fields(self) -> tuple[str, ...]:
        return self.columns
 

    def create(
        self,
        schema: Schema,
        **kwargs,
    ) -> h5py.Dataset:
        """Create a 1D-datasate"""
        return super().create(schema=schema, shape=(0,), **kwargs)


    def open(
        self,
        schema: Schema | None = None,
        *,
        strict: bool = True,
    ) -> h5py.Dataset:
        """Open the underlying HDF5 dataset and populate lookup caches.

        Parameters
        ----------
        schema : Schema or None, optional
            Schema used to validate schema identity. If None, uses the cached schema from the last successful create/open.
        strict : bool, optional
            If True, enforce strict schema identity checks.

        Returns
        -------
        h5py.Dataset
            Open dataset handle.

        Raises
        ------
        ValueError
            If `mode` is invalid.
        """

        ds = super().open(schema, strict=strict)
        return ds

    def add(
        self,
        rows: list[dict],
        *,
        defaults: dict[str, object] | None = None,
    ) -> None:
        """Append rows to the registry (append-only insert).

        Parameters
        ----------
        rows : list[dict]
            Records to insert. Keys matching dataset fields will be written.
   
        """
        pass 
      

    def update(self, rows: list[dict], *, warn_missing=True) -> None:
        """Update the metadata of existing registered entities.

        rows : list[dict]
            Records to update. Keys matching dataset fields will be updated.
        warn_missing : bool
            Whether or not to warn user if a provided ``id`` is not in the registry.

        Returns
        -------
        None

        """
        pass
 




    def to_dataframe(self, arr=None, **_):
        # If caller didn't pass an array, slice from disk.
        if arr is None:
            ds = self._ds()
            # Empty dataset: return empty DataFrame with on-disk columns.
            if ds.shape[0] == 0:
                return pd.DataFrame(columns=list(ds.dtype.names or ()))
            arr = ds[:]

        # If caller passed an empty structured array, preserve column names.
        if (
            arr is not None
            and getattr(getattr(arr, "dtype", None), "names", None) is not None
            and getattr(arr, "size", 0) == 0
        ):
            return pd.DataFrame(columns=list(arr.dtype.names or ()))

        arr = decode_arr(arr)

        # `decode_arr` may return the array unchanged for empty structured arrays.
        if (
            isinstance(arr, np.ndarray)
            and getattr(arr.dtype, "names", None) is not None
            and arr.size == 0
        ):
            return pd.DataFrame(columns=list(arr.dtype.names or ()))

        df = pd.DataFrame(arr)

        # decode any bytes that survived into object columns
        for col in df.columns:
            if df[col].dtype == object:
                df[col] = df[col].map(
                    lambda x: (
                        x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else x
                    )
                )

        return df

    def to_csv(self,
        fname: str,
        *,
        delimiter=",",
        **kwargs):
        """Save a table to a csv file.
        """
        df = self.to_dataframe(**kwargs)
        df.to_csv(fname, sep=delimiter, **kwargs)

    def validate(self, *, strict: bool = True) -> None:
        """Validate basic invariants of the provenance dataset.

        Currently checks only:
        - dataset exists
        - dataset is 1D resizable
        - dtype is a compound dtype (recommended)

        This stays lightweight by design.
        """
        if not self.exists():
            if strict:
                raise RuntimeError(f"Dataset does not exist at {self.path}")
            return

        ds = self._ds()
        if ds.ndim != 1:
            raise ValueError(f"Provenance dataset must be 1D; found ndim={ds.ndim}")
        if (
            ds.maxshape is not None
            and ds.maxshape[0] is not None
            and ds.maxshape[0] < ds.shape[0]
        ):
            raise ValueError(
                "Provenance dataset maxshape is inconsistent with current shape"
            )

        # Not strictly required, but provenance should almost always be a compound dtype.
        if ds.dtype.names is None:
            msg = f"Provenance dataset dtype is not compound at {self.path}"
            if strict:
                raise ValueError(msg)



    def get(
        self, selector: Callable | None=None
    ) -> pd.DataFrame:
        """Retrieve  rows as a pandas DataFrame.

        Parameters
        ----------
        selector : object, optional
            If None, return all rows under `mode`. If callable, it is treated as a
            predicate and applied to a materialized DataFrame. Otherwise it is
            interpreted according to `by`.


        Returns
        -------
        pandas.DataFrame
            Selected rows.


        """


        df = self.to_dataframe()

        if selector is None:
            return df

        
        return df[selector(df)]






    

    