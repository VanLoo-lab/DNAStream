"""dnastream.registry

This module implements an HDF5-backed append-only registry with label/id resolution.

The primary entry point is :class:`~dnastream.registry.Registry`.
"""

from __future__ import annotations
from typing import Mapping
import uuid
import getpass
from datetime import datetime, timezone
import numpy as np
import h5py
import warnings
from ._h5base import H5Dataset
from .utils import norm_key, as_str, as_str_vec, _qualname, decode_arr

import pandas as pd
import collections.abc as cabc
from .schema import Schema
from .constants import (
    VALID_MODES,
    BY_OPTIONS,
    REGISTRY_SPINE,
    LABEL_SCALARS,
    ID_SCALARS,
    EVENTS,
    SCOPE,
)


class Registry(H5Dataset):
    """HDF5-backed append-only registry with active/inactive history.

    A `Registry` stores rows in an HDF5 compound dataset. Each row has a unique
    identifier (`id`) and may have a `label`. The `active` flag indicates which
    row is the current active representative for a label.

    The registry is append-only: inserts add new rows; updates are represented by
    changing state flags (e.g., deactivating an older row).

    Parameters
    ----------
    parent : h5py.Group
        Parent HDF5 group. Must have name ``"/registry"``.
    name : str
        Dataset name within ``/registry``.

    Raises
    ------
    ValueError
        If `parent` is not the ``/registry`` group.

    Notes
    -----
    The registry maintains in-memory caches to accelerate resolution from labels
    and ids to row indices. Cache contents depend on the view `mode` used when
    opening or resolving.
    """

    def __init__(self, parent: h5py.Group, name: str):
        """Create a new `Registry` wrapper.

        Parameters
        ----------
        parent : h5py.Group
            Parent group. Must be ``"/registry"``.
        name : str
            Dataset name.

        Notes
        -----
        This constructor does not create the dataset. It only binds to an HDF5 parent
        and initializes internal caches.
        """
        if parent.name != "/registry":
            raise ValueError("Registries must live within the registry group.")
        super().__init__(parent, name)

        # initialize caches
        self._cache_valid = False
        self._label_to_id = None
        self._id_to_label = None
        self._label_to_idx = None
        self._id_to_idx = None

    def __len__(self):
        return self._ds().shape[0]

    def __repr__(self) -> str:
        n = "?"  # avoid disk access if dataset isn't available
        try:
            n = len(self)
        except Exception:
            pass
        return f"Registry(name={self.name!r}, path={self.path!r}, " f"n={n})"

    def __iter__(self):
        for row in self._ds():
            yield decode_arr(row)

    def __contains__(self, item) -> bool:
        """
        Returns true if uuid id in registry

        Notes
        -----
        This checks membership within the entire registy whether active/inactive
        """
        if item is None:
            return False

        is_scalar, items = self._validate_id_selector(item)
        if not is_scalar:
            raise ValueError(f"Item {item} must be a scalar value.")

        item = norm_key(as_str(items[0]))

        self._load_cache()
        return item in self._id_to_idx

    def __getitem__(self, item) -> dict[str, object]:
        """Return a registry row (as a dict) by id.

        Parameters
        ----------
        item
            Registry id (UUID, str, or bytes). Must be a scalar.

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
        """
        if item is None:
            raise KeyError("None is not a valid registry id.")

        is_scalar, items = self._validate_id_selector(item)
        if not is_scalar:
            raise TypeError(
                f"Registry id lookup expects a scalar; got {type(item).__name__}."
            )

        rid = norm_key(as_str(items[0]))

        self._load_cache()
        return decode_arr(self._row_by_id(rid))

    @property
    def columns(self) -> tuple[str, ...]:
        ds = self.open()
        return tuple(ds.dtype.names or ())

    @property
    def fields(self) -> tuple[str, ...]:
        cols = self.columns
        spine = set(REGISTRY_SPINE)
        return tuple(c for c in cols if c not in spine)

    def create(
        self,
        schema: Schema,
        **kwargs,
    ) -> h5py.Dataset:
        """Create a registry dataset and log"""
        super().create(schema=schema, shape=(0,), **kwargs)
        self._emit(
            EVENTS.CREATE,
            fn=_qualname(self.create),
            schema_version=schema.version,
            schema_hash=schema.hash(),
        )

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
        self._load_cache()
        return ds

    def add(
        self,
        rows: list[dict],
        activate_newest: bool = True,
        allow_duplicate_labels=False,
    ) -> None:
        """Append rows to the registry (append-only insert).

        Parameters
        ----------
        rows : list[dict]
            Records to insert. Keys matching dataset fields will be written.
        activate_newest : bool, optional
            Collision policy when an incoming label matches an *active* existing label.

            - If True (default): deactivate the existing active row and keep the new row
              active.
            - If False: keep the existing row active and insert the new row as inactive.
        allow_duplicate_labels : bool, optional
            Whether to allow duplicate labels to be added to the registry

            - If True: data associated wtih labels that collide with existing labels will be added
            to the registry and activation will be set according to `activate_newest`.
            - If False (default), data associated with existing labels will not be added to the registry.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If label computation fails, produced labels are empty, or duplicates exist
            within the incoming batch.

        Notes
        -----
        The following fields are auto-populated:

        - ``id``: UUID4 string
        - ``idx``: on-disk row index
        - ``active``: activatation status of entity
        - ``created_at``: UTC ISO 8601 timestamp ending in ``"Z"``
        - ``created_by``: current username
        """
        if len(rows) == 0:
            return

        ds = self.open()
        n0 = ds.shape[0]
        names = ds.dtype.names

        n_add = len(rows)

        # Compute labels (or None) and find collisions against *active* existing rows
        if not self.schema.label_required:
            labels = None
        else:
            labels = as_str_vec([self.schema.make_label(row) for row in rows])
            if any(not lab for lab in labels):
                raise ValueError("Empty label produced for at least one input row.")

        collisions: list[tuple[int, int]] = []
        if labels is not None:
            collisions = self._detect_label_collisions(labels, active_only=True)

        if not allow_duplicate_labels:
            duplicate_label_idx = [x for x, _ in collisions]
            n_add = len(rows) - len(duplicate_label_idx)
        else:
            duplicate_label_idx = []

        if n_add == 0:
            return

        protected = set(REGISTRY_SPINE)

        now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        user = getpass.getuser()

        # Resize and allocate block

        block = np.empty(n_add, dtype=ds.dtype)

        for name in ds.dtype.names:
            kind = ds.dtype[name].kind
            if kind in ("O", "U"):  # vlen/unicode strings
                block[name] = ""
            elif kind == "S":  # fixed-width bytes
                block[name] = b""
            elif kind == "b":  # bool
                block[name] = False
            elif kind in ("i", "u"):  # ints
                block[name] = 0
            elif kind == "f":  # floats
                block[name] = np.nan
            else:
                # as a last resort; but try to avoid leaving garbage
                block[name] = None
        # Defaults

        # Fill in registry spine iaw the registry contract
        block["id"] = [str(uuid.uuid4()) for _ in range(n_add)]

        if labels is None:
            block["label"] = ""
        else:
            block["label"] = labels

        block["active"] = True
        block["created_at"] = now
        block["created_by"] = user
        block["modified_at"] = now
        block["modified_by"] = user
        block["idx"] = np.arange(n0, n0 + n_add)

        # Apply collision policy
        if collisions and allow_duplicate_labels:
            if activate_newest:
                # Deactivate existing active rows
                self.deactivate_labels([labels[new_i] for new_i, _ in collisions])
            else:
                # Keep existing active: mark new colliding rows inactive+deprecated
                new_idxs = [new_i for new_i, _old_i in collisions]
                block["active"][new_idxs] = False

        # Fill from input
        for i, row in enumerate(rows):
            # skip duplicates
            if i in duplicate_label_idx:
                continue
            for name in names:
                if name in protected:
                    continue

                if name in row and row[name] is not None:
                    field = self.schema.field(name)
                    value = row[name]
                    if field.validator is not None:
                        field.validator(value)

                    if ds.dtype[name].kind in ("O", "U"):
                        block[name][i] = str(value)
                    elif ds.dtype[name].kind == "S":
                        block[name][i] = str(value).encode("utf-8")
                    else:
                        block[name][i] = value

        ds.resize((n0 + n_add,))
        ds[n0 : n0 + n_add] = block

        self._invalidate_cache()
        self._emit(EVENTS.APPEND, _qualname(self.add), n_add=n_add)

    def update(self, rows: list[dict], *, warn_missing=True) -> None:
        """Update the metadata of existing registered entities.

        rows : list[dict]
            Records to update. Keys matching dataset fields will be updated.
        warn_missing : bool
            Whether or not to warn user if a provided ``id`` is not in the registry.

        Returns
        -------
        None

        Notes
        ------
        Due to Registry policy allowing duplicate labels, the ``id`` must be provided as part of the
        dictionary updates. use `resolve_ids(...)` to look up the ids for each activated label. Use `find_ids(...)`
        to look up ``id`` in the event of duplicate labels.


        Raises
        -------
        ValueError : if ``id`` is not a key in a dictionary in one of the input rows to be edited.
        """
        ds = self.open()
        names = set(ds.dtype.names or ())

        if not rows:
            return

        now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        user = getpass.getuser()

        non_editable = set(REGISTRY_SPINE)
        if self.schema.label_from is not None:
            non_editable |= set(self.schema.label_from)

        self._load_cache()

        records_to_update = []
        indices = []
        for row in rows:
            if "id" not in row:
                raise ValueError(
                    "'id' field not found. Must be present in dictionary to update registry!"
                )

            rid = norm_key(as_str(row["id"]))

            idx = self._idx_for_id(rid)

            if idx is None:

                if warn_missing:

                    self._warn_missing([rid])

                continue

            rec = self._row(idx)

            changed = False
            for key, val in row.items():
                if key == "id":
                    continue
                if key in non_editable:
                    warnings.warn(
                        f"Field {key} is non-editable, metadata for {rid!r} and {key} not modified."
                    )
                    continue
                if key not in names:
                    continue

                # Encode / coerce to match HDF5 dtype expectations
                field = self.schema.field(key)
                if field.validator is not None:
                    field.validator(val)
                kind = ds.dtype[key].kind
                if kind == "S":  # fixed-width bytes
                    rec[key] = ("" if val is None else str(val)).encode("utf-8")
                elif kind in ("U", "O"):  # unicode/object
                    rec[key] = "" if val is None else str(val)
                else:
                    rec[key] = val

                changed = True

            if changed:

                rec["modified_by"] = (
                    str(user).encode("utf-8")
                    if ds.dtype["modified_by"].kind == "S"
                    else str(user)
                )
                rec["modified_at"] = (
                    str(now).encode("utf-8")
                    if ds.dtype["modified_at"].kind == "S"
                    else str(now)
                )

                indices.append(idx)
                records_to_update.append(rec)

        indices = np.asarray(indices, np.int64)
        records = np.asarray(records_to_update, dtype=ds.dtype)

        sorted_idx = np.argsort(np.asarray(indices, np.int64))
        sorted_indices = indices[sorted_idx]
        sorted_records = records[sorted_idx]

        ds[sorted_indices] = sorted_records

        self._emit(
            EVENTS.MODIFY, _qualname(self.update), n_modified=len(sorted_records)
        )

    def resolve_ids(self, labels, *, missing=np.nan) -> np.ndarray:
        """Resolve label(s) to ACTIVE id(s) '``.

        Parameters
        ----------
        labels
            Label or iterable of labels.
        missing : object, optional
            Value returned for missing labels.

        Returns
        -------
        numpy.ndarray
            Array of ids aligned to input order.

        Notes
        ------
        Only resolves ids for active labels. Use `find_ids(...)` to find all all ids associated with a given label.

        """
        # validate the input

        is_scalar, sel = self._validate_label_selector(labels)

        self._load_cache()
        return self._resolve_from_map(
            sel,
            self._label_to_id,
            missing=missing,
            key_norm=norm_key,
            is_scalar=is_scalar,
        )

    def resolve_labels(self, ids, *, missing=np.nan) -> np.ndarray:
        """Resolve id(s) to label(s).

        Parameters
        ----------
        ids
            Id or iterable of ids.
        missing : object, optional
            Value returned for missing ids.

        Returns
        -------
        numpy.ndarray
            Array of labels aligned to input order.

        """

        is_scalar, sel = self._validate_id_selector(ids)

        self._load_cache()
        return self._resolve_from_map(
            sel,
            self._id_to_label,
            missing=missing,
            key_norm=norm_key,
            is_scalar=is_scalar,
        )

    def find_ids(self, labels, *, mode: str = "all") -> dict[str, np.ndarray]:
        """Return all ids matching each label under the requested mode.

        Parameters
        ----------
        labels
            Label or iterable of labels.
        mode : {"all", "non_active"}, optional
            Which subset of rows to consider.

        Returns
        -------
        dict[str, numpy.ndarray]
            Mapping from normalized label to an array of matching ids.

        Notes
        -----
        Unlike :meth:`resolve_ids`, this can return multiple ids per label.
        """
        return self._find_multiple(labels, mode=mode, return_field="id")

    def activate_labels(self, labels, activate_newest=True, warn_missing=True):
        """Activate registry entries by label.

        Parameters
        ----------
        labels
            Label or iterable of labels.
        activate_newest : bool, optional
            When no active row exists for a label, choose the newest historical row if
            True else the oldest.
        warn_missing : bool, optional
            If True, emit warnings for missing labels.
        """
        _, labels = self._validate_label_selector(labels)

        self._activate(
            labels,
            by="label",
            activate_newest=activate_newest,
            warn_missing=warn_missing,
        )

        self._emit(EVENTS.DESIGNATE, _qualname(self.activate_labels))

    def activate_ids(self, ids, *, warn_missing=True):
        """Activate registry entries by id.

        Parameters
        ----------
        ids
            Id or iterable of ids.
        warn_missing : bool, optional
            If True, emit warnings for missing ids.
        """
        _, ids = self._validate_id_selector(ids)
        self._activate(ids, by="id", warn_missing=warn_missing)

        self._emit(EVENTS.DESIGNATE, _qualname(self.activate_ids))

    def deactivate_labels(self, labels, *, warn_missing=True) -> None:
        """Deactivate registry entries by label.

        Parameters
        ----------
        labels
            Label or iterable of labels.
        warn_missing : bool, optional
            If True, emit warnings for missing labels.

        Notes
        -----
        Deactivates all entities that match the corresponding label, even in cases with multiple entries. For specific,
        entities, it is recommended to use `deactivate_ids()`

        """
        _, labels = self._validate_label_selector(labels)
        self._deactivate(labels, by="label", warn_missing=warn_missing)
        self._emit(EVENTS.DESIGNATE, _qualname(self.deactivate_labels))

    def deactivate_ids(self, ids, *, warn_missing=True) -> None:
        """Deactivate registry entries by id.

        Parameters
        ----------
        ids
            Id or iterable of ids.
        warn_missing : bool, optional
            If True, emit warnings for missing ids.
        """
        _, ids = self._validate_id_selector(ids)
        self._deactivate(ids, by="id", warn_missing=warn_missing)
        self._emit(EVENTS.DESIGNATE, _qualname(self.deactivate_ids))

    def to_dataframe(
        self,
        arr: np.ndarray | None = None,
        *,
        mode: str = "all",
        **kwargs,
    ) -> pd.DataFrame:
        """Materialize registry rows as a pandas DataFrame.

        Parameters
        ----------
        arr : numpy.ndarray or None, optional
            If provided, convert this already-materialized structured array to a
            DataFrame (no additional disk reads). This exists so callers like
            `get()` can slice on-disk first and then convert.
        mode : {"active_only", "non_active", "all"}, optional
            Which subset of rows to include when `arr` is None.

        Returns
        -------
        pandas.DataFrame
            DataFrame view of the requested rows.

        Notes
        -----
        This reads from disk when `arr` is None.
        """
        # Fast path: caller already sliced an array.
        if arr is not None:
            return super().to_dataframe(arr)

        mode = self._validate_mode(mode)
        ds = self._ds()

        if mode == "all":
            arr = ds[:]
        else:
            active = ds["active"][:].astype(bool)
            mask = active if mode == "active_only" else ~active
            idx = np.nonzero(mask)[0].astype(np.int64)
            arr = ds[idx] if idx.size else np.zeros((0,), dtype=ds.dtype)

        return super().to_dataframe(arr)

    def validate(self, *, strict: bool = True) -> None:
        """Validate registry invariants.

        Parameters
        ----------
        strict : bool, optional
            If True (default), raise :class:`ValueError` on any violation. If False,
            emit warnings instead.

        Returns
        -------
        None

        Raises
        ------
        ValueError
            If `strict=True` and any invariants are violated.

        Notes
        -----
        The following invariants are checked:

        - Dataset path is ``/registry/<name>``.
        - Required spine fields are present.
        - Ids are unique and valid UUIDs.
        - ``idx`` equals the on-disk row index and is in range.
        - Active labels are unique.
        """
        ds = self._ds()
        errors: list[str] = []

        ### check the path is properly formatted
        path = ds.name
        if path != f"/registry/{self.name}":
            errors.append(
                f"Registry path: {path} is invalid. Should be /registry/<name>"
            )

        # check required fields are in dtype names
        names = set(ds.dtype.names or ())

        # registry spine are required fields for any registyr

        missing_fields = set(REGISTRY_SPINE) - names
        if missing_fields:
            errors.append(f"Missing required field(s): {sorted(missing_fields)}")

        ##### check ids #####
        if "id" not in missing_fields:
            ids = ds["id"][:]
            norm_ids = [norm_key(x) for x in ids]
            # check uniqueness
            if len(set(norm_ids)) != len(norm_ids):
                errors.append("Duplicate ids found in registry.")

            # check valid uuids
            for i, rid in enumerate(norm_ids):

                try:
                    uuid.UUID(rid)
                except Exception:
                    errors.append(f"Row {i} has invalid UUID {rid!r}")
        else:
            errors.append(
                "id validation checks not performed due to missing field 'id'"
            )

        #### check row idx####
        if "idx" not in missing_fields:
            n = ds.shape[0]
            for i, idx in enumerate(ds["idx"][:]):
                try:
                    idx_i = int(idx)
                except Exception:
                    errors.append(f"Row {i} has non-integer idx {idx!r}")
                    continue
                if idx_i != i:
                    errors.append(f"Row {i} does not match idx {idx_i}")
                if not (0 <= idx_i < n):
                    errors.append(f"Row {i} idx out of range: {idx_i} (n={n})")
        else:
            errors.append(
                "idx validation checks not performed due to missing field 'idx'"
            )

        #### check labels#####
        if "active" not in missing_fields and "label" not in missing_fields:
            active = ds["active"][:].astype(bool)
            labels = ds["label"][:]
            seen = {}
            for i, (lab, ok) in enumerate(zip(labels, active)):
                if not ok:
                    continue
                key = norm_key(lab)
                if key in seen:
                    errors.append(
                        f"Duplicate active label '{key}' at rows {seen[key]} and {i}"
                    )
                seen[key] = i
        else:
            errors.append(
                "label validation checks not performed due to missing field 'label' or 'active'"
            )

        if errors:
            msg = (
                f"Registry '{self.name}' failed validation "
                f"with {len(errors)} error(s):\n"
                + "\n".join(f"  - {e}" for e in errors)
            )
            if strict:
                raise ValueError(msg)
            else:
                warnings.warn(msg, stacklevel=2)

    def get(
        self, selector=None, *, by: str | None = None, mode: str = "active_only"
    ) -> pd.DataFrame:
        """Retrieve registry rows as a pandas DataFrame.

        Parameters
        ----------
        selector : object, optional
            If None, return all rows under `mode`. If callable, it is treated as a
            predicate and applied to a materialized DataFrame. Otherwise it is
            interpreted according to `by`.
        by : {"label", "id", "idx"} or None, optional
            How to interpret `selector` when it is not callable.
        mode : {"active_only", "non_active", "all"}, optional
            Registry view to query.

        Returns
        -------
        pandas.DataFrame
            Selected rows.

        Raises
        ------
        ValueError
            If `by` is invalid or missing when required.

        Notes
        -----
        For `by != "idx"`, this resolves selectors to on-disk indices using caches
        and may scan the dataset in history modes.
        """
        mode = self._validate_mode(mode)

        if selector is None:
            return self.to_dataframe(mode=mode)

        if callable(selector):
            # this path inherently needs a DataFrame of *something* to call the predicate on
            # you can either: (a) keep it as "loads df", or (b) disallow callables for lazy mode
            df = self.to_dataframe(mode=mode)
            return df[selector(df)]

        if by is None or by not in {"id", "label", "idx"}:
            raise ValueError(
                "Must pass by='label', by='id', or by='idx' when selector is not callable."
            )

        if by == "idx":
            idx = np.asarray(selector, dtype=object)
            if idx.ndim == 0:
                idx = np.array([idx], dtype=object)
        else:
            if by == "label":
                is_scalar, sel = self._validate_label_selector(selector)
            else:
                is_scalar, sel = self._validate_id_selector(selector)

            idx = self._select_indices(selector=sel, by=by, mode=mode)

        # normalize scalar -> 1D
        idx = np.asarray(idx, dtype=object)

        # drop missing and convert to int
        idx = idx[~pd.isna(idx)].astype(int)

        # IMPORTANT: h5py fancy indexing is happiest with sorted unique indices.
        # Also avoids duplicated disk reads.
        order = np.argsort(idx)
        idx_sorted = idx[order]

        ds = self._ds()
        arr = ds[idx_sorted]  # <-- on-disk slice only
        df = self.to_dataframe(arr=arr)

        # restore original requested order (optional but usually expected)
        inv = np.empty_like(order)
        inv[order] = np.arange(order.size)
        return df.iloc[inv].reset_index(drop=True)

    def _invalidate_cache(self) -> None:
        """Invalidate all cached lookup mappings.

        Notes
        -----
        Caches are rebuilt on demand by :meth:`open` / :meth:`_load_cache`.
        """
        self._cache_valid = False
        self._label_to_id = None
        self._id_to_label = None
        self._label_to_idx = None
        self._id_to_idx = None

    @staticmethod
    def _validate_mode(mode: str) -> str:
        """Validate and normalize registry view modes.

        Parameters
        ----------
        mode : str or None
            Desired cache/view mode. If None, defaults to ``"all"``.

        Returns
        -------
        str
            Normalized lowercase mode.

        Raises
        ------
        ValueError
            If `mode` is not one of ``{"active_only", "non_active", "all"}``.
        """
        if mode is None:
            mode = "all"

        mode = str(mode).strip().lower()

        if mode not in VALID_MODES:
            raise ValueError(
                f"Invalid mode {mode!r}. Expected one of {sorted(VALID_MODES)}."
            )
        return mode

    @staticmethod
    def _validate_by(by: str) -> str:
        """Validate and normalize selector interpretation.

        Parameters
        ----------
        by : str
            Selector interpretation.

        Returns
        -------
        str
            Normalized lowercase selector type.

        Raises
        ------
        ValueError
            If `by` is None or not one of ``{"label", "id"}``.
        """
        if by is None:
            raise ValueError(f"by cannot be None; expected one of {sorted(BY_OPTIONS)}")
        by = str(by).strip().lower()
        if by not in BY_OPTIONS:
            raise ValueError(
                f"Invalid by={by!r}. Expected one of {sorted(BY_OPTIONS)}."
            )
        return by

    def _load_cache(self) -> None:
        """Build in-memory caches for fast id/label resolution.



        Raises
        ------
        ValueError
            If required fields are missing or if duplicates violate registry invariants.

        Notes
        -----
        - In ``"active_only"`` mode, label resolution assumes labels are unique among
          active rows.
        - In other modes, label-based caches are not built because labels may map to
          multiple rows.
        """

        if self._cache_valid:
            return
        ds = self._ds()
        names = ds.dtype.names

        required = set(REGISTRY_SPINE)
        missing = required - set(names)
        if missing:
            raise ValueError(
                f"Registry '{self.name}' missing fields: {sorted(missing)}"
            )

        ids = ds["id"][:]

        labels = ds["label"][:]
        idxs = ds["idx"][:]
        active = ds["active"][:].astype(bool)

        label_to_id = {}
        id_to_label = {}
        id_to_idx = {}
        label_to_idx = {}

        for lab, rid, idx, ok in zip(labels, ids, idxs, active):

            labn = norm_key(lab)
            ridn = norm_key(rid)

            if ridn in id_to_idx:
                raise ValueError(f"Duplicate id '{ridn}' in registry '{self.name}'")

            id_to_idx[ridn] = int(idx)
            id_to_label[ridn] = labn

            if ok:
                if labn in label_to_idx:
                    raise ValueError(
                        f"Duplicate active label '{labn}' in registry '{self.name}'"
                    )
                label_to_idx[labn] = int(idx)
                label_to_id[labn] = ridn

        self._label_to_id = label_to_id
        self._label_to_idx = label_to_idx
        self._id_to_idx = id_to_idx
        self._id_to_label = id_to_label
        self._cache_valid = True

    def _set_state(
        self,
        indices,
        *,
        active: bool | None = None,
    ) -> None:
        """Set registry state flags for the given row indices.

        Parameters
        ----------
        indices : int or list[int] or numpy.ndarray
            Row indices to update.
        active : bool or None, optional
            If provided, set the `active` field to this value for all selected rows.

        Returns
        -------
        None

        Notes
        -----
        This mutates the on-disk dataset fields in-place.
        """
        ds = self._ds()

        if active is None:
            return

        # h5py: make scalar index work with fancy indexing
        _, idx = self._normalize_selector(indices, scalar_types=(int, np.integer))

        if len(idx) == 0:
            return

        idx = np.sort(np.asarray(idx, np.int64))

        rows = ds[idx]
        rows["active"] = bool(active)
        ds[idx] = rows

        after = ds["active"][idx]
        if not np.all(after == active):
            raise RuntimeError(f"Correct state '{active}' not set on disk")

    def _detect_label_collisions(
        self,
        labels: list[str] | None,
        *,
        active_only: bool = True,
    ) -> list[tuple[int, int]]:
        if labels is None:
            return []

        # Normalize first; duplicates should be defined on normalized keys
        labels_norm = [norm_key(as_str(lab)) for lab in labels]
        if len(set(labels_norm)) != len(labels_norm):
            raise ValueError(
                "Duplicate labels supplied in input data (after normalization). "
                f"Check registry '{self.name}' schema for label definitions."
            )

        ds = self._ds()
        if ds.shape[0] == 0:
            return []

        names = ds.dtype.names
        if "label" not in names:
            raise ValueError(f"Registry '{self.name}' missing required field 'label'.")

        existing_labels = ds["label"][:]

        if active_only:
            # spine guarantees active exists; fail loudly if not
            if "active" not in names:
                raise ValueError(
                    f"Registry '{self.name}' missing required field 'active'."
                )
            active_mask = ds["active"][:].astype(bool)
        else:
            active_mask = None

        # Build label -> existing row index mapping
        # - active_only=True: duplicates are a registry invariant violation
        # - active_only=False: duplicates are history; keep newest row index
        existing_index: dict[str, int] = {}
        for i, lab in enumerate(existing_labels):
            if active_mask is not None and not bool(active_mask[i]):
                continue
            key = norm_key(as_str(lab))

            if active_only and key in existing_index:
                raise ValueError(
                    f"Registry '{self.name}' contains multiple active rows with label '{key}'."
                )

            existing_index[key] = i  # newest wins (esp. for history mode)

        collisions: list[tuple[int, int]] = []
        for new_i, key in enumerate(labels_norm):
            old_i = existing_index.get(key)
            if old_i is not None:
                collisions.append((new_i, old_i))

        return collisions

    def _row(self, idx: int | np.integer) -> np.void:
        """Return the on-disk row at integer index `idx`.

        Notes
        -----
        This is an internal helper. It validates bounds against the dataset shape.
        """
        ds = self._ds()

        if idx is None:
            raise TypeError("idx cannot be None")

        try:
            idx_i = int(idx)
        except (TypeError, ValueError) as e:
            raise TypeError(f"idx must be an integer; got {type(idx).__name__}.") from e

        if idx_i < 0 or idx_i >= len(self):
            raise IndexError(
                f"Row index {idx_i} out of bounds for dataset of length {len(self)}"
            )

        return ds[idx_i]

    def _idx_for_id(self, id) -> int | None:
        """Return row index for a given id, or None if missing."""
        # Ensure id->idx cache exists (id caches are valid for any mode).
        self._load_cache()
        key = norm_key(as_str(id))
        return self._id_to_idx.get(key, None)

    def _row_by_id(self, id) -> np.void:
        """Return the on-disk row for a given id.

        Raises
        ------
        KeyError
            If the id is not present.
        """
        idx = self._idx_for_id(id)
        if idx is None:
            raise KeyError(f"id {id!r} not found in Registry '{self.name}'.")
        return self._row(idx)

    def _resolve_from_map(
        self,
        keys,
        mapping: Mapping[str, object] | None,
        *,
        missing=np.nan,
        key_norm=None,
        is_scalar=False,
    ):
        """Resolve one key or many keys using a mapping.

        Parameters
        ----------
        keys
            Scalar key or iterable of keys.
        mapping : Mapping[str, object] or None
            Mapping used to resolve each key. If None, all lookups return `missing`.
        missing : object, optional
            Value used when a key is not found.
        key_norm : callable or None, optional
            Optional normalization function applied to each key before lookup.
        scalar_types : tuple[type, ...], optional
            Types treated as scalar (non-iterable) keys.

        Returns
        -------
        object or numpy.ndarray
            If `keys` is scalar, returns a scalar. Otherwise returns a 1D object array.

        Notes
        -----
        This is a pure in-memory operation and does not touch HDF5.
        """

        if mapping is None:
            out = np.asarray([missing] * len(keys), dtype=object)
            return out[0] if is_scalar else out

        if key_norm is None:
            it = (mapping.get(k, missing) for k in keys)
        else:
            it = (mapping.get(key_norm(k), missing) for k in keys)

        out = np.fromiter(it, dtype=object, count=len(keys))
        return out[0] if is_scalar else out

    def _select_indices(
        self,
        selector,
        by: str,
        *,
        mode: str = "active_only",
        missing=np.nan,
        warn_missing: bool = True,
        allow_many: bool = True,
        select_newest: bool = True,
    ) -> np.ndarray:
        """Resolve a user selector into row indices.

        Responsibilities
        ----------------
        - Normalize selector (scalar vs iterable)
        - Apply view mode filtering (active_only / non_active / all)
        - Warn for missing selections (optional)
        - Return indices sorted/unique for h5py fancy indexing

        Parameters
        ----------
        selector
            Scalar or iterable of ids/labels/indices depending on `by`.
        by : {'id','label'}
            How to interpret `selector`.
        mode : {'active_only','non_active','all'}
            Registry view to apply.
        missing
            Sentinel used for unresolved lookups.
        warn_missing
            If True, emit warnings for unresolved requested selections.
        allow_many
            If True and by='label' with mode != 'active_only', return *all* matching
            indices for each label (history). If False, pick a single winner per label.
        select_newest
            When allow_many=False and there is no active row, choose newest (max idx)
            if True else oldest (min idx).

        Returns
        -------
        np.ndarray
            Sorted unique integer row indices.
        """
        mode = self._validate_mode(mode)
        by = self._validate_by(by)

        # Ensure caches exist for this view. Note: id->idx is valid in any mode.
        self._load_cache()
        ds = self._ds()

        def _empty() -> np.ndarray:
            return np.asarray([], dtype=np.int64)

        def _unique_int(x) -> np.ndarray:
            if x is None:
                return _empty()
            arr = np.asarray(x, dtype=np.int64)
            if arr.size == 0:
                return _empty()
            return np.unique(arr)

        # ---- by == 'id' -------------------------------------------------
        if by == "id":
            _, sel = self._validate_id_selector(selector)
            idx = self._resolve_from_map(
                sel,
                self._id_to_idx,
                missing=missing,
                key_norm=norm_key,
            )
            idx = np.asarray(idx, dtype=object)
            miss_mask = pd.isna(idx)

            if warn_missing and miss_mask.size:
                self._warn_missing(
                    [req for req, m in zip(sel, miss_mask.tolist()) if m]
                )

            all_idx = idx[~miss_mask].astype(np.int64)
            if all_idx.size == 0:
                return _empty()

            if mode == "all":
                return _unique_int(all_idx)

            active_state = ds["active"][all_idx].astype(bool)
            if mode == "active_only":
                return _unique_int(all_idx[active_state])
            else:  # non_active
                return _unique_int(all_idx[~active_state])

        # ---- by == 'label' ----------------------------------------------
        # Fast path: active_only has unique label->idx cache
        if mode == "active_only":
            _, sel = self._validate_label_selector(selector)
            idx = self._resolve_from_map(
                sel,
                self._label_to_idx,
                missing=missing,
                key_norm=norm_key,
            )
            idx = np.asarray(idx, dtype=object)
            miss_mask = pd.isna(idx)

            if warn_missing and miss_mask.size:
                self._warn_missing(
                    [req for req, m in zip(sel, miss_mask.tolist()) if m]
                )

            out = idx[~miss_mask].astype(np.int64)
            return _unique_int(out)

        # History modes: labels may map to multiple rows.
        _, sel = self._validate_label_selector(selector)
        want = [norm_key(x) for x in sel]
        hits = self._find_multiple(want, mode=mode, return_field="idx")

        missing_labels = [
            lab for lab, arr in hits.items() if arr is None or len(arr) == 0
        ]
        if warn_missing:
            self._warn_missing(sorted(missing_labels))

        if allow_many:
            flat: list[int] = []
            for arr in hits.values():
                if arr is None:
                    continue
                flat.extend([int(i) for i in np.asarray(arr, dtype=object).tolist()])
            return _unique_int(flat)

        # allow_many=False: choose a single winner per label
        winners: list[int] = []
        for lab_norm, idxs_arr in hits.items():
            if idxs_arr is None:
                continue
            idxs = [int(i) for i in np.asarray(idxs_arr, dtype=object).tolist()]
            if len(idxs) == 0:
                continue

            # Prefer an already-active row if present (should be unique)
            active_state = ds["active"][idxs].astype(bool)
            actives = [i for i, ok in zip(idxs, active_state) if ok]

            if len(actives) > 1:
                raise ValueError(
                    f"Duplicate active labels found for label {lab_norm!r} in Registry '{self.name}'."
                )

            if len(actives) == 1:
                winners.append(int(actives[0]))
            else:
                winners.append(int(max(idxs) if select_newest else min(idxs)))

        return _unique_int(winners)

    def _find_multiple(
        self,
        labels,
        *,
        return_field: str = "id",
        mode: str = "all",
    ) -> dict[str, np.ndarray]:
        """Return all matches for each requested label under the requested mode.

        Parameters
        ----------
        labels
            Label or iterable of labels.
        return_field : {"id", "idx"}, optional
            Which field to return for each matching row.
        mode : {"all", "non_active"}, optional
            Which subset of rows to consider.

        Returns
        -------
        dict[str, numpy.ndarray]
            Mapping from normalized label to an array of matching values.

        Raises
        ------
        ValueError
            If `mode` or `return_field` are invalid, or required fields are missing.

        Notes
        -----
        This method scans the dataset once and can return multiple values per label.
        """
        if mode not in {"all", "non_active"}:
            raise ValueError("mode must be 'all' or 'non_active'")

        if return_field not in {"id", "idx"}:
            raise ValueError("return_field must be 'id' or 'idx'")

        ds = self._ds()
        names = ds.dtype.names
        if "label" not in names:
            raise ValueError(f"Registry '{self.name}' must contain 'label'")
        if return_field not in names:
            raise ValueError(f"Registry '{self.name}' must contain '{return_field}'")

        # normalize input
        scalar_types = (str, bytes, np.bytes_, np.str_)
        is_scalar = isinstance(labels, scalar_types)
        label_list = [labels] if is_scalar else list(labels)
        want = [norm_key(x) for x in label_list]

        # pre-seed outputs (so every requested label has an entry)
        out_lists: dict[str, list[object]] = {w: [] for w in want}
        if ds.shape[0] == 0 or not out_lists:
            return {k: np.asarray([], dtype=object) for k in out_lists}

        # row mask
        if mode == "non_active" and "active" in names:
            mask = ~ds["active"][:].astype(bool)
        else:
            mask = np.ones(ds.shape[0], dtype=bool)

        labs = ds["label"][:]
        vals = ds[return_field][:]

        want_set = set(out_lists.keys())

        # single scan through registry
        if return_field == "id":
            for lab, val, ok in zip(labs, vals, mask):
                if not ok:
                    continue
                key = norm_key(lab)
                if key in want_set:
                    out_lists[key].append(norm_key(val))
        else:  # idx
            for lab, val, ok in zip(labs, vals, mask):
                if not ok:
                    continue
                key = norm_key(lab)
                if key in want_set:
                    out_lists[key].append(int(val))

        return {k: np.asarray(v, dtype=object) for k, v in out_lists.items()}

    def _warn_missing(self, missing_items) -> None:

        if len(missing_items) == 0:
            return
        warnings.warn(
            f"Requested selection(s) {missing_items!r} not found in Registry '{self.name}'.",
            stacklevel=2,
        )

    @staticmethod
    def _normalize_selector(
        selector, scalar_types=(str, bytes, np.bytes_, np.str_)
    ) -> tuple[bool, list]:
        """Normalize a selector to a list.

        Returns (is_scalar, sel_list). Treat numpy 0-d arrays as scalars (unwrap).
        Reject mappings (dict) because iterating them yields keys (usually a bug).
        """
        # unwrap numpy 0-d arrays
        if isinstance(selector, np.ndarray) and selector.shape == ():
            selector = selector.item()

        # reject mappings explicitly
        if isinstance(selector, cabc.Mapping):
            raise ValueError(
                "Selector must be a scalar or an iterable of scalars; mappings (e.g., dict) are not supported."
            )

        # scalar path
        if isinstance(selector, scalar_types):
            return True, [selector]

        # iterable path
        try:
            return False, list(selector)
        except TypeError as e:
            raise ValueError(
                f"Selector must be a scalar or an iterable of scalars; got {type(selector).__name__}."
            ) from e

    @staticmethod
    def _validate_selector_scalar_types(selector, *, allowed_scalar_types, what: str):
        # allow numpy 0-d arrays but unwrap them
        if isinstance(selector, np.ndarray) and selector.shape == ():
            selector = selector.item()

        # scalar path
        if isinstance(selector, allowed_scalar_types):
            return

        # if it’s not one of the allowed scalar types, it might be an iterable;
        # let _normalize_selector handle it later. BUT reject scalar-like non-iterables.
        # (Also reject common “scalar but not allowed” immediately.)
        if np.isscalar(selector):
            raise ValueError(
                f"{what} must be one of {allowed_scalar_types}; got {type(selector).__name__}."
            )

    @staticmethod
    def _validate_selector_elements(seq, *, allowed_scalar_types, what: str):
        for x in seq:
            if isinstance(x, np.ndarray) and x.shape == ():
                x = x.item()
            if not isinstance(x, allowed_scalar_types):
                raise ValueError(
                    f"Each {what} must be one of {allowed_scalar_types}; got {type(x).__name__}."
                )

    def _validate_id_selector(self, selector) -> tuple[bool, list]:
        self._validate_selector_scalar_types(
            selector, allowed_scalar_types=ID_SCALARS, what="id"
        )
        is_scalar, sel = self._normalize_selector(selector, scalar_types=ID_SCALARS)
        self._validate_selector_elements(
            sel, allowed_scalar_types=ID_SCALARS, what="id"
        )

        return is_scalar, sel

    def _validate_label_selector(self, selector) -> tuple[bool, list]:
        self._validate_selector_scalar_types(
            selector, allowed_scalar_types=LABEL_SCALARS, what="label"
        )
        is_scalar, sel = self._normalize_selector(selector, scalar_types=LABEL_SCALARS)
        self._validate_selector_elements(
            sel, allowed_scalar_types=LABEL_SCALARS, what="label"
        )

        return is_scalar, sel

    def _activate(
        self, selector, by, *, activate_newest: bool = True, warn_missing=True
    ) -> None:
        """Activate registry entries by label or id.

        This enforces the invariant that at most one row per normalized label is
        active at any time.

        Parameters
        ----------
        selector
            A single label/id or an iterable of labels/ids.
        by : {'label','id'}
            Select whether `selector` is interpreted as labels or ids.
        activate_newest : bool
            When a label has no active row, choose which historical row to
            activate. If True, activates the most recently appended row
            (largest row index). If False, activates the oldest row.

        Notes
        -----
        - If a label already has an active row, that row is kept active and all
          other rows with the same label are deactivated.
        - Missing requested labels/ids emit warnings.
        - This mutates the on-disk `active` field and invalidates caches.
        """
        by = self._validate_by(by)

        if by == "label":
            _, selector = self._validate_label_selector(selector)
        else:
            _, selector = self._validate_id_selector(selector)

        activations = self._select_indices(
            selector,
            by=by,
            mode="all",
            allow_many=False,
            select_newest=activate_newest,
            warn_missing=warn_missing,
        )

        if by == "id":
            selector = self.resolve_labels(selector)

            selector = selector.tolist()
        all_indices = self._select_indices(
            selector, by="label", mode="all", allow_many=True, warn_missing=False
        )

        deactivations = sorted(set(all_indices) - set(activations))

        # Apply activated state changes
        if activations.size:
            self._set_state(indices=activations, active=True)
        if deactivations:
            self._set_state(indices=deactivations, active=False)

        self._invalidate_cache()

    def _deactivate(self, selector, by, *, warn_missing=True) -> None:
        """Deactivate registry entries by label or id.

        Parameters
        ----------
        selector
            A single label/id or an iterable of labels/ids.
        by : {'label','id'}
            Select whether `selector` is interpreted as labels or ids.

        Notes
        -----
        - If there are multiple registry entities associated with a label, all are deactivated.
        - Missing requested labels/ids emit warnings.
        - This mutates the on-disk `active` field and invalidates caches.
        """
        by = self._validate_by(by)
        indices = self._select_indices(
            selector, by=by, mode="all", warn_missing=warn_missing
        )

        if indices.size:

            self._set_state(indices=indices, active=False)

        self._invalidate_cache()

    def _emit(self, event: str, fn: str, **payload):
        return super()._emit(SCOPE.REGISTRY, event, fn, **payload)
