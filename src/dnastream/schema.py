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

from dataclasses import dataclass, field
from typing import Callable, Any
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
    dtype: Any  # e.g. np.dtype("i4") or "U50"
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
        Callable ``label_builder(parts) -> str`` where ``parts`` is a
        ``tuple[str, ...]`` extracted from ``label_from``. If None, defaults to
        joining parts with ``"|"``.
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
    label_builder: Callable[[tuple[str, ...]], str] | None = None
    label_normalizer: Callable[[str], str] | None = None
    _field_by_name: dict[str, Field] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Compute derived attributes (dtype, field index) for a frozen dataclass."""
        fb = {f.name: f for f in self.fields}

        if len(fb) != len(self.fields):
            raise ValueError("Duplicate field names in schema.")

        if self.label_from is None and self.label_required:
            raise ValueError("'label_from' cannot be 'None' when label is required.")

        if self.label_from and self.label_required:
            for x in self.label_from:
                if str(x) not in fb:
                    raise ValueError(f" label_from Field {x!r} not found in Schema.")

        object.__setattr__(self, "dtype", np.dtype(list(self.get_spec())))
        object.__setattr__(self, "_field_by_name", fb)

    def field(self, name: str) -> Field:
        """Return the :class:`~dnastream.schema.Field` for a given name.

        Parameters
        ----------
        name : str
            Field name.

        Returns
        -------
        Field
            The corresponding Field object.

        Raises
        ------
        KeyError
            If ``name`` is not present in the schema.
        """
        try:
            return self._field_by_name[name]
        except KeyError:
            raise KeyError(f"Field {name!r} not in schema.") from None

    def required_names(self) -> tuple[str, ...]:
        """Return the names of required fields.

        Returns
        -------
        tuple[str, ...]
            Names of fields with ``required=True`` in schema order.
        """
        return tuple(f.name for f in self.fields if f.required)

    def make_label(self, row: dict[str, Any]) -> str:
        """Build a label string from an input record.

        The label pipeline is:

        1. Extract values from ``row`` using ``label_from``.
        2. Convert extracted values to strings.
        3. Call ``label_builder(parts)`` (or default join with ``"|"``).
        4. Apply ``label_normalizer`` if provided.

        Parameters
        ----------
        row : dict[str, Any]
            Input record.

        Returns
        -------
        str
            The constructed label.

        Raises
        ------
        ValueError
            If ``label_from`` is None.
        KeyError
            If a key in ``label_from`` is missing from ``row``.
        """
        if self.label_from is None:
            raise ValueError("Schema.label_from is None; cannot build label")

        parts = tuple(str(row[k]) for k in self.label_from)

        if self.label_builder is None:
            label = "|".join(parts)
        else:
            if len(parts) > 1:
                # Support both: label_builder(*parts) and label_builder(parts)

                label = self.label_builder(*parts)
            else:

                label = self.label_builder(parts)

        if self.label_normalizer is not None:
            label = self.label_normalizer(label)

        return label

    def get_spec(self):
        """Return the dtype specification tuple for the schema.

        Returns
        -------
        tuple[tuple[str, Any], ...]
            Tuple of ``(name, dtype)`` pairs in schema order.
        """
        return tuple((field.name, field.dtype) for field in self.fields)

    def get_names(self):
        """Return field names in schema order.

        Returns
        -------
        tuple[str, ...]
            Field names.
        """
        return tuple(field.name for field in self.fields)

    def schema_pairs(self):
        """Return a JSON-friendly list of ``(name, dtype_str)`` pairs.

        Returns
        -------
        list[tuple[str, str]]
            Pairs suitable for serialization.
        """
        return [(k, str(v)) for k, v in self.get_spec()]

    def get_payload(self):
        """Return the canonical JSON payload used for hashing.

        Returns
        -------
        str
            Compact JSON string encoding of ``schema_pairs``.
        """
        return json.dumps(
            sorted([(k, str(v)) for k, v in self.get_spec()], key=lambda x: x[0]),
            separators=(",", ":"),
        )

    def json_pairs(self):
        """Return the canonical JSON string for this schema.

        Returns
        -------
        str
            JSON string encoding of the schema pairs.
        """
        return self.get_payload()

    def hash(self):
        """Compute a stable SHA-256 hash of the schema payload.

        Returns
        -------
        str
            Hex digest of SHA-256 over the UTF-8 encoded payload.
        """
        payload = self.get_payload().encode("utf-8")
        return hashlib.sha256(payload).hexdigest()
