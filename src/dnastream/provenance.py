from __future__ import annotations
import json
from dataclasses import dataclass
from datetime import datetime, timezone
import getpass
import uuid
from typing import Any, Callable, Iterable

import h5py
import numpy as np

from ._h5base import H5Dataset
from .constants import EVENTS, Event, Scope

from .utils import _qualname


class Provenance(H5Dataset):
    """HDF5-backed append-only provenance log.

    A `Provenance` stores rows in an HDF5 compound dataset. Each row has a unique
    identifier (`id`). The log is append-only.

    Notes
    -----
    This wrapper assumes the underlying dataset already exists (created via
    :meth:`H5Dataset.create`) and that the schema for this dataset is cached in
    ``self._schema`` by the surrounding `DNAStream` machinery.
    """

    def __init__(self, parent: h5py.Group, name: str):
        if parent.name != "/provenance":
            raise ValueError("Logs must live within the provenance group.")
        super().__init__(parent, name)

    def __len__(self) -> int:
        return self._ds().shape[0]

    def __repr__(self) -> str:
        n = "?"  # avoid disk access if dataset isn't available
        try:
            n = len(self)
        except Exception:
            pass
        return f"Provenance(name={self.name!r}, path={self.path!r}, n={n})"

    def __iter__(self) -> Iterable[np.void]:
        for row in self._ds():
            yield row

    @staticmethod
    def _init_scalar(dtype: np.dtype) -> np.ndarray:
        """Create a scalar structured array and initialize string-like fields."""
        row = np.zeros((), dtype=dtype)
        if dtype.names is not None:
            for name in dtype.names:
                kind = dtype[name].kind
                if kind in ("O", "U"):
                    row[name] = ""
                elif kind == "S":
                    row[name] = b""
        return row

    def add(
        self,
        scope: Scope,
        event: Event,
        dataset: str,
        fn: str,
        **payload: Any,
    ) -> None:
        """Append a provenance record.

        Parameters
        ----------
        event
            One of the allowed provenance event strings (see ``dnastream.constants.EVENTS``).
        dataset
            Full HDF5 path of the dataset being acted on (e.g. ``"/registry/variant"``).
        fn
            Optional callable responsible for the action; recorded as ``source``.
        **payload
            Arbitrary event metadata; stored as a JSON string in ``info``.
        """

        if event not in EVENTS:
            raise ValueError(
                f"Invalid event {event!r}. Expected one of {sorted(EVENTS)}"
            )

        timestamp = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
        user = getpass.getuser()
        entry_id = str(uuid.uuid4())
        source = fn

        # Keep this robust: provenance should never fail just because a payload value
        # isn't JSON-serializable.
        info = json.dumps(payload, sort_keys=True, default=str) if payload else ""

        ds = self.open()
        n0 = len(self)
        ds.resize((n0 + 1,))

        dt = ds.dtype

        # We expect a compound dtype with named fields per PROVENANCE_LOG_SCHEMA.
        if dt.names is None:
            raise TypeError(
                f"Provenance dataset at {self.path} must have a compound dtype with named fields"
            )

        row = self._init_scalar(dt)

        if "id" in dt.names:
            row["id"] = entry_id
        if "timestamp" in dt.names:
            row["timestamp"] = timestamp
        if "user" in dt.names:
            row["user"] = user
        if "scope" in dt.names:
            row["scope"] = scope
        if "event" in dt.names:
            row["event"] = event
        if "dataset" in dt.names:
            row["dataset"] = dataset
        if "source" in dt.names:
            row["source"] = source
        if "info" in dt.names:
            row["info"] = info

        ds[n0] = row

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

    def register_hook(self, *args, **kwargs):
        raise RuntimeError("Hooks are not supported on provenance datasets.")

    def _emit(self, *args, **kwargs):
        raise RuntimeError(
            "Provenance must not emit hooks/events (would recurse). "
            "If you meant to log something, call ds.provenance('changes').add(...) "
            "from the *caller*, not from inside Provenance."
        )
