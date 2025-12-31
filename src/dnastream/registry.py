from __future__ import annotations
from typing import Mapping
import uuid
import getpass
from datetime import datetime, timezone
import numpy as np
import h5py
import warnings
from ._h5base import H5Dataset
from .utils import norm
import pandas as pd

_VALID_MODES = {"active_only", "all", "non_active"}


class Registry(H5Dataset):
    def __init__(self, parent: h5py.Group, name: str):
        if parent.name != "/registry":
            raise ValueError("Registries must live within the registry group.")
        super().__init__(parent, name)

        # initialize caches
        self._cache_valid = False
        self._cache_mode = None  # None / "active_only" / "non_active" / "all"
        self._label_to_id = None
        self._id_to_label = None
        self._label_to_idx = None
        self._id_to_idx = None

    def _invalidate_cache(self) -> None:
        self._cache_valid = False
        self._cache_mode = None
        self._label_to_id = None
        self._id_to_label = None
        self._label_to_idx = None
        self._id_to_idx = None

    def open(
        self, expected_schema: dict | None = None, *, strict: bool = True
    ) -> h5py.Dataset:
        """
        Opens a handle to the registry and loads the caches
        """
        ds = super.open(expected_schema, strict=strict)
        self._load_cache(active_only=False)
        return ds

    def _load_cache(self, *, cache_mode="active_only"):
        if self._cache_valid and self._cache_mode == cache_mode:
            return
        ds = self._ds()
        names = ds.dtype.names

        required = {"id", "label", "idx", "active"}
        missing = required - set(names)
        if missing:
            raise ValueError(
                f"Registry '{self.name}' missing fields: {sorted(missing)}"
            )

        ids = ds["id"][:]
        labels = ds["label"][:]
        idxs = ds["idx"][:]
        active = ds["active"][:].astype(bool)

        if cache_mode == "active_only":
            mask = active
        elif cache_mode == "non_active":
            mask = np.logical_not(active)
        else:
            mask = np.ones(ds.shape[0], dtype=bool)

        label_to_id = {}
        id_to_label = {}
        id_to_idx = {}
        label_to_idx = {}

        for lab, rid, idx, ok in zip(labels, ids, idxs, mask):
            if not ok:
                continue
            labn = norm(lab)
            ridn = norm(rid)
            if labn in label_to_idx:
                if cache_mode == "active_only":
                    raise ValueError(
                        f"Duplicate active label '{labn}' in registry '{self.name}'"
                    )
            if ridn in id_to_idx:
                raise ValueError(f"Duplicate id '{ridn}' in registry '{self.name}'")
            label_to_idx[labn] = int(idx)
            id_to_idx[ridn] = int(idx)
            label_to_id[labn] = ridn
            id_to_label[ridn] = labn

        # don't populate potential many to one caches (used find_*() instead)
        if cache_mode != "active_only":
            self._label_to_idx = None
            self.label_to_id = None
        else:
            self._label_to_id = label_to_id
            self._label_to_idx = label_to_idx

        self._id_to_idx = id_to_idx
        self._id_to_label = id_to_label
        self._cache_valid = True
        self._cache_mode = cache_mode

    def _set_state(
        self,
        indices,
        *,
        active: bool | None = None,
    ) -> None:
        """
        Set registry state flags for the given row indices.

        Parameters
        ----------
        indices : int | list[int] | np.ndarray
            Row indices to update.
        active : bool | None
            If provided, set `active` to this value.

        """
        ds = self._ds()
        names = ds.dtype.names

        if active is None:
            return

        # h5py: make scalar index work with fancy indexing
        if isinstance(indices, (int, np.integer)):
            idx = [int(indices)]
        else:
            idx = indices

        if active is not None:
            if "active" not in names:
                raise ValueError(f"Registry '{self._name}' has no 'active' field")
            ds["active"][idx] = bool(active)

    def _label_normalizer(self, rows: list[dict], schema: dict) -> list[str] | None:
        label_required = schema.get("label_required", False)

        if not label_required:
            return None

        labels: list[str] = []
        label_from = schema.get("label_from", None)
        label_normalizer_fn = schema.get("label_normalizer", None)

        if not callable(label_normalizer_fn):
            raise ValueError(
                "schema.label_normalizer must be callable when label_from is set."
            )

        for idx, r in enumerate(rows):

            if "label" in r:
                warnings.warn(f"label {r['label']} overwritten")

            try:
                args = [r[f] for f in label_from]
            except KeyError as e:
                raise ValueError(
                    f"Row {idx} missing required label field {e!s}"
                ) from None

            lab = label_normalizer_fn(*args)

            if not lab:
                raise ValueError(f"Row {idx} produced an empty label.")

            labels.append(lab)

        return labels

    def _detect_label_collisons(
        self,
        labels: list[str] | None,
        *,
        active_only: bool = True,
    ) -> list[tuple[int, int]]:
        """

        Detect label collisions between incoming rows and existing registry entries.

        Parameters
        ----------
        labels : list[str] | None
            Proposed labels for the incoming rows. If None, returns [].
        active_only : bool
            If True (default), only consider existing rows with active==True
            when checking for collisions.

        Returns
        -------
        list[tuple[int, int]]
            List of (new_row_idx, existing_row_idx) pairs for each incoming label
            that collides with an existing label under the chosen filter.

        Raises
        ------
        ValueError
            - If duplicate labels are present within the incoming batch.
            - If the registry contains multiple *active* rows with the same label.

        Notes
        -----
        This method does not mutate registry state.
        Collision resolution is the responsibility of the caller.
        """
        if labels is None:
            return []

        # 1) duplicates within incoming batch are ambiguous raise Value Error
        if len(set(labels)) != len(labels):
            raise ValueError(
                f"Duplicate labels supplied in input data. "
                f"Check registry '{self._name}' schema for label definitions."
            )

        ds = self._ds()
        if ds.shape[0] == 0:
            return []

        names = ds.dtype.names
        have_active = "active" in names

        existing_labels = ds["label"][:]
        if active_only and have_active:
            active_mask = ds["active"][:].astype(bool)
        else:
            active_mask = np.ones(ds.shape[0], dtype=bool)

        # 2) Build label -> existing index mapping for considered rows
        existing_index: dict[str, int] = {}
        for i, (lab, is_active) in enumerate(zip(existing_labels, active_mask)):
            if not is_active:
                continue
            s = lab.decode("utf-8") if isinstance(lab, (bytes, np.bytes_)) else str(lab)

            if s in existing_index:
                # two active rows with same label = data integrity bug
                raise ValueError(
                    f"Registry '{self._name}' contains multiple active rows with label '{s}'."
                )
            existing_index[s] = i

        # 3) Return collision pairs
        collisions: list[tuple[int, int]] = []
        for new_i, lab in enumerate(labels):
            key = str(lab)
            old_i = existing_index.get(key)
            if old_i is not None:
                collisions.append((new_i, old_i))

        return collisions

    def add(self, rows: list[dict], schema: dict, keep_new: bool = True) -> None:
        """Append-only insert.

        Parameters
        ----------
        rows : list[dict]
            Records to insert.
        schema : dict
            Compiled schema from `dnastream.schemas.compile_schema`.
        keep_new : bool
            Collision policy when an incoming label matches an *active* existing label.
            If True (default): deactivate,the existing active row; keep the new row active.
            If False: keep existing active; insert the new row as inactive.

        """
        ds = self.open(schema)
        n0 = ds.shape[0]
        names = ds.dtype.names

        n_add = len(rows)
        if n_add == 0:
            return []

        now = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        user = getpass.getuser()

        # Compute labels (or None) and find collisions against *active* existing rows
        labels = self._label_normalizer(rows, schema)
        collisions: list[tuple[int, int]] = []
        if labels is not None:
            collisions = self._detect_label_collisons(labels, active_only=True)

        # Resize and allocate block
        ds.resize((n0 + n_add,))
        block = np.zeros(n_add, dtype=ds.dtype)

        # Defaults
        if "id" in names:
            block["id"] = [str(uuid.uuid4()) for _ in range(n_add)]
        if "label" in names:
            if labels is None:
                block["label"] = ""
            else:
                block["label"] = np.asarray(labels, dtype=object)

        if "active" in names:
            block["active"] = True
        if "created_at" in names:
            block["created_at"] = now
        if "created_by" in names:
            block["created_by"] = user

        if "idx" in names:
            block["idx"] = np.arange(n0, n0 + n_add)

        # Apply collision policy
        if collisions:
            if keep_new:
                # Deactivate existing active rows
                old_idxs = [old_i for _new_i, old_i in collisions]
                self._set_state(old_idxs, active=False)
            else:
                # Keep existing active: mark new colliding rows inactive+deprecated
                new_idxs = [new_i for new_i, _old_i in collisions]
                if "active" in names:
                    block["active"][new_idxs] = False

        # Fill from input
        for i, row in enumerate(rows):
            for name in names:
                if name in row and row[name] is not None:
                    if ds.dtype[name].kind in ("O", "U"):
                        block[name][i] = str(row[name])
                    elif ds.dtype[name].kind == "S":
                        block[name][i] = str(row[name]).encode("utf-8")
                    else:
                        block[name][i] = row[name]

        ds[n0 : n0 + n_add] = block

        self._invalidate_cache()

    def _resolve_from_map(
        self,
        keys,
        mapping: Mapping[str, object] | None,
        *,
        missing=np.nan,
        key_norm=None,
        scalar_types=(str, bytes, np.bytes_, np.str_),
    ):
        """
        Resolve one key or many keys using mapping.get(key, missing).

        Notes
        -----
        - Does NOT touch HDF5. Assumes mapping represents the desired view (active_only, etc.).
        - If mapping is None, returns all missing.
        """
        # handle numpy 0-d arrays as scalar
        if isinstance(keys, np.ndarray) and keys.shape == ():
            keys = keys.item()

        is_scalar = isinstance(keys, scalar_types)
        key_list = [keys] if is_scalar else list(keys)

        if mapping is None:
            out = np.asarray([missing] * len(key_list), dtype=object)
            return out[0] if is_scalar else out

        if key_norm is None:
            it = (mapping.get(k, missing) for k in key_list)
        else:
            it = (mapping.get(key_norm(k), missing) for k in key_list)

        out = np.fromiter(it, dtype=object, count=len(key_list))
        return out[0] if is_scalar else out

    def resolve_id(self, labels, *, mode="active_only", missing=np.nan) -> np.ndarray:
        if mode != "active_only":
            raise ValueError(
                "Only mode='active_only' is allowed when resolving ids from labels. "
                "Use find_id(...) for non-active lookups."
            )
        self._load_cache(cache_mode=mode)
        return self._resolve_from_map(
            labels,
            self._label_to_id,
            missing=missing,
            key_norm=norm,
        )

    def resolve_label(self, ids, *, mode="active_only", missing=np.nan) -> np.ndarray:
        self._load_cache(cache_mode=mode)
        return self._resolve_from_map(
            ids,
            self._id_to_label,
            missing=missing,
            key_norm=norm,
        )

    def resolve_idx_from_ids(
        self, ids, *, mode="active_only", missing=np.nan
    ) -> np.ndarray:
        self._load_cache(cache_mode=mode)
        return self._resolve_from_map(
            ids,
            self._id_to_idx,
            missing=missing,
            key_norm=norm,
        )

    def resolve_idx_from_labels(
        self, labels, *, mode="active_only", missing=np.nan
    ) -> np.ndarray:
        if mode != "active_only":
            raise ValueError(
                f"Only mode active_only allowed when resolving indices from labels."
            )
        self._load_cache(cache_mode=mode)
        return self._resolve_from_map(
            labels,
            self._label_to_idx,
            missing=missing,
            key_norm=norm,
        )

    def _find_multiple(
        self,
        labels,
        *,
        return_field: str = "id",
        mode: str = "all",
    ) -> dict[str, np.ndarray]:
        """
        Return all ids/idx matching each label under the requested mode.
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
        want = [norm(x) for x in label_list]

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
                key = norm(lab)
                if key in want_set:
                    out_lists[key].append(norm(val))
        else:  # idx
            for lab, val, ok in zip(labs, vals, mask):
                if not ok:
                    continue
                key = norm(lab)
                if key in want_set:
                    out_lists[key].append(int(val))

        return {k: np.asarray(v, dtype=object) for k, v in out_lists.items()}

    def find_ids(self, labels, *, mode: str = "all") -> dict[str, np.ndarray]:
        """
        Return all ids matching each label under the requested mode.

        Notes
        -----
        Unlike resolve_id(), this can return multiple ids per label when labels are
        not unique (e.g., inactive history).
        """
        return self._find_multiple(labels, mode=mode, return_field="id")

    def find_rows(self, labels, *, mode: str = "all") -> dict[str, np.ndarray]:
        """
        Return all row indices matching each label under the requested mode.

        Notes
        -----
        Unlike resolve_id(), this can return multiple ids per label when labels are
        not unique (e.g., inactive history).
        """
        return self._find_multiple(labels, mode=mode, return_field="idx")

    def activate(self, labels, mode="non_active", activate_newest=True):
        """
        Given a list of labels, activate sets the activate state to True
        if the requested label is found. If collison


        Notes
        -----
        Warns users if requested labels to activate are not found in the registry.



        """
        collisions = self._detect_label_collisons(labels, mode="active_only")

        ids_to_activate = self.resolve_idx(labels, mode=mode, missing="")

        self._set_state()

    @staticmethod
    def _validate_mode(mode: str) -> str:
        """Validate and normalize registry view modes."""
        if mode is None:
            mode = "all"

        mode = str(mode).strip().lower()

        if mode not in _VALID_MODES:
            raise ValueError(
                f"Invalid mode {mode!r}. Expected one of {sorted(_VALID_MODES)}."
            )
        return mode

    def to_dataframe(self, mode="all") -> pd.DataFrame:
        mode = self._validate_mode(mode)
        ds = self._ds()

        if mode == "all":
            arr = ds[:]
        else:

            active = ds["active"][:].astype(bool)
            mask = active if mode == "active_only" else ~active
            idx = np.nonzero(mask)[0].astype(np.int64)
            arr = ds[idx] if idx.size else np.zeros((0,), dtype=ds.dtype)

        return self._to_dataframe(arr)

    def validate(self, *, strict: bool = True) -> None:
        """
        Validate registry invariants.

        Raises
        ------
        ValueError
            If strict=True and any violations are found.
        Warns
        -----
        UserWarning
            If strict=False and violations are found.

        Notes
        -----
        - checks that ids are unique
        - checks that registry path is valid
        - checks that fields include required fields id, label, idx and active
        - checks row idx are valid, i.e., between 0...size -1 and idx equals row index
        - checks ids are non-null and valid uuids
        - checks labels are unique in active only
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
        required = {"active", "id", "idx", "label"}
        missing_fields = required - names
        if missing_fields:
            errors.append(f"Missing required field(s): {sorted(missing_fields)}")

        ##### check ids #####
        if "id" not in missing_fields:
            ids = ds["id"][:]
            norm_ids = [norm(x) for x in ids]
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
                key = norm(lab)
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
