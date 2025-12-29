import uuid
import getpass
from datetime import datetime, timezone
import numpy as np
import h5py
import json
import warnings
from ._h5base import H5Dataset


class Registry(H5Dataset):
    def __init__(self, handle: h5py.File, grp: str, name: str):
        super().__init__(handle, grp, name)
        self.id_col = "id"

    def deprecate(self, labels=None, deactivated_only=True):
        """
        Udpates the registry by deprecating all labels or a specified list of labels.


        Parameters
        ----------
        labels : list[str]
            labels to deprecate

        deactivated_only: boolean
            if true, only labels that are currently deactivated are deprecated


        """

        pass

    def _set_state(
        self,
        indices,
        *,
        active: bool | None = None,
        deprecated: bool | None = None,
    ) -> None:
        """
        Set registry state flags for the given row indices.

        Parameters
        ----------
        indices : int | list[int] | np.ndarray
            Row indices to update.
        active : bool | None
            If provided, set `active` to this value.
        deprecated : bool | None
            If provided, set `deprecated` to this value.

        Invariants
        ----------
        - If active is set True and deprecated is not explicitly provided,
        deprecated is forced False (can't be both active and deprecated).
        """
        ds = self._ds()
        names = ds.dtype.names

        if active is None and deprecated is None:
            return

        # h5py: make scalar index work with fancy indexing
        if isinstance(indices, (int, np.integer)):
            idx = [int(indices)]
        else:
            idx = indices

        if active is True and deprecated is None:
            deprecated = False

        if active is not None:
            if "active" not in names:
                raise ValueError(f"Registry '{self._name}' has no 'active' field")
            ds["active"][idx] = bool(active)

        if deprecated is not None:
            if "deprecated" not in names:
                raise ValueError(f"Registry '{self._name}' has no 'deprecated' field")
            ds["deprecated"][idx] = bool(deprecated)

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

    def _label_validator(
        self,
        labels: list[str] | None,
        *,
        active_only: bool = True,
    ) -> list[tuple[int, int]]:
        """Find collisions between incoming labels and existing registry labels.

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

    def add(self, rows: list[dict], schema: dict, keep_new: bool = True):
        """Append-only insert.

        Parameters
        ----------
        rows : list[dict]
            Records to insert.
        schema : dict
            Compiled schema from `dnastream.schemas.compile_schema`.
        keep_new : bool
            Collision policy when an incoming label matches an *active* existing label.
            If True (default): deactivate+deprecate the existing active row; keep the new row active.
            If False: keep existing active; insert the new row as inactive+deprecated.

        Returns
        -------
        list[tuple[int, int]]
            Collision pairs (new_row_idx, existing_row_idx).
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
            collisions = self._label_validator(labels, active_only=True)

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
        if "deprecated" in names:
            block["deprecated"] = False
        if "created_at" in names:
            block["created_at"] = now
        if "created_by" in names:
            block["created_by"] = user

        # Apply collision policy
        if collisions:
            if keep_new:
                # Deactivate+deprecate existing active rows
                old_idxs = [old_i for _new_i, old_i in collisions]
                self._set_state(old_idxs, active=False, deprecated=True)
            else:
                # Keep existing active: mark new colliding rows inactive+deprecated
                new_idxs = [new_i for new_i, _old_i in collisions]
                if "active" in names:
                    block["active"][new_idxs] = False
                if "deprecated" in names:
                    block["deprecated"][new_idxs] = True

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
        return collisions
