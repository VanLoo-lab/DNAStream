"""dnastream.schema

Lightweight schema objects used to define and validate HDF5-backed datasets.

This module defines two immutable dataclasses:

- :class:`~dnastream.schema.Field`: a single column specification.
- :class:`~dnastream.schema.Schema`: an ordered collection of fields plus
  schema metadata (version, hashing helpers, optional label construction).

Notes
-----
These classes are intentionally small and dependency-light. They are designed
to be serializable *by identity* (e.g., version + hash + json pairs) while
keeping non-serializable callables (validators/label functions) in-memory.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Any, Optional
import json
import hashlib
import numpy as np


@dataclass(frozen=True)
class Field:
    """A single column specification within a :class:`~dnastream.schema.Schema`.

    Parameters
    ----------
    name : str
        Field/column name.
    dtype : Any
        NumPy/h5py dtype specification for the column, e.g. ``np.int64`` or
        ``h5py.string_dtype("utf-8")``.
    required : bool, default True
        If True, the field must be present (and non-None) in incoming records
        during validation.
    validator : callable or None, optional
        Optional callable ``validator(value) -> None`` that should raise an
        exception if ``value`` is invalid.

    Notes
    -----
    Validation semantics (e.g., how ``required`` is enforced) are implemented
    by the consumer (e.g., a Registry) rather than by :class:`Field` itself.
    """

    name: str
    dtype: Any
    required: bool = True
    validator: Callable[[Any], None] | None = None  # raise on failure


@dataclass(frozen=True)
class Schema:
    """An immutable registry schema.

    A :class:`Schema` defines the ordered set of fields for a registry and
    provides helpers for dtype construction, stable identity (hashing), and
    optional label generation.

    Parameters
    ----------
    fields : tuple[Field, ...]
        Ordered field specifications.
    version : str
        Human-readable schema version string (e.g., ``"1.0.0"``).
    label_from : tuple[str, ...] or None, optional
        If provided, label components are extracted from these keys in an input
        record and passed to ``label_builder``.
    label_required : bool, default False
        If True, a non-empty label is required for each inserted record.
    label_builder : callable or None, optional
        Callable used to construct a label from the extracted components.

        By default, labels are built by joining components with ``"|"``.

        If provided, the builder is called as ``label_builder(*parts)`` where
        ``parts`` are the stringified values extracted from ``label_from``.
        For backward compatibility, builders that accept a single tuple argument
        may also be used (they will be called as ``label_builder(parts)``).
    label_normalizer : callable or None, optional
        Callable ``label_normalizer(label) -> str`` applied after
        ``label_builder``.

    Attributes
    ----------
    dtype : numpy.dtype
        Structured dtype derived from ``fields``.

    Notes
    -----
    ``Schema`` is frozen (immutable). Derived attributes (e.g., ``dtype`` and
    the name->Field index) are computed in ``__post_init__``.

    Callables (validators and label functions) are intentionally not serialized
    into HDF5. Consumers should persist schema identity (e.g., version + hash)
    and reconstruct the runtime Schema in code.
    """

    fields: tuple[Field, ...]
    version: str
    dtype: np.dtype = field(init=False, repr=False)

    label_from: tuple[str, ...] | None = None
    label_required: bool = False
    label_builder: Optional[Callable[..., str]] = None
    label_normalizer: Callable[[str], str] | None = None

    _field_by_name: dict[str, Field] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        fb = {f.name: f for f in self.fields}

        if len(fb) != len(self.fields):
            raise ValueError("Duplicate field names in schema.")

        if self.label_from is None and self.label_required:
            raise ValueError("'label_from' cannot be 'None' when label is required.")

        if self.label_from and self.label_required:
            for x in self.label_from:
                if str(x) not in fb:
                    raise ValueError(f"label_from field {x!r} not found in schema.")

        object.__setattr__(self, "dtype", np.dtype(list(self.get_spec())))
        object.__setattr__(self, "_field_by_name", fb)

    def field(self, name: str) -> Field:
        try:
            return self._field_by_name[name]
        except KeyError:
            raise KeyError(f"Field {name!r} not in schema.") from None

    def required_names(self) -> tuple[str, ...]:
        return tuple(f.name for f in self.fields if f.required)

    # -----------------------------
    # Docs/helpers
    # -----------------------------
    def field_summary(self) -> list[dict[str, Any]]:
        """Return a JSON-friendly summary of fields for documentation."""
        out: list[dict[str, Any]] = []
        for f in self.fields:
            v = f.validator
            out.append(
                {
                    "name": f.name,
                    "dtype": str(f.dtype),
                    "required": bool(f.required),
                    "validator": (
                        getattr(v, "__name__", None) if v is not None else None
                    ),
                }
            )
        return out

    def to_markdown_table(self) -> str:
        """Render fields as a Markdown table: name | dtype | required | validator."""
        rows = self.field_summary()
        cols = ("name", "dtype", "required", "validator")
        header = "| " + " | ".join(cols) + " |"
        sep = "| " + " | ".join(["---"] * len(cols)) + " |"
        lines = [header, sep]
        for r in rows:
            lines.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
        return "\n".join(lines)

    # -----------------------------
    # Labeling
    # -----------------------------
    def make_label(self, row: dict[str, Any]) -> str:
        if self.label_from is None:
            raise ValueError("Schema.label_from is None; cannot build label")

        parts = tuple(str(row[k]) for k in self.label_from)

        if self.label_builder is None:
            label = "|".join(parts)
        else:
            try:
                label = self.label_builder(*parts)
            except TypeError:
                label = self.label_builder(parts)

        if self.label_normalizer is not None:
            label = self.label_normalizer(label)

        return label

    # -----------------------------
    # Identity / dtype helpers
    # -----------------------------
    def get_spec(self) -> tuple[tuple[str, Any], ...]:
        return tuple((field.name, field.dtype) for field in self.fields)

    def get_names(self) -> tuple[str, ...]:
        return tuple(field.name for field in self.fields)

    def schema_pairs(self) -> list[tuple[str, str]]:
        return [(k, str(v)) for k, v in self.get_spec()]

    def get_payload(self) -> str:
        return json.dumps(
            sorted([(k, str(v)) for k, v in self.get_spec()], key=lambda x: x[0]),
            separators=(",", ":"),
        )

    def json_pairs(self) -> str:
        return self.get_payload()

    def hash(self) -> str:
        payload = self.get_payload().encode("utf-8")
        return hashlib.sha256(payload).hexdigest()
