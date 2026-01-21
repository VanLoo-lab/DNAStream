from dataclasses import dataclass, field
from typing import Callable, Any
import json
import hashlib
import numpy as np


@dataclass(frozen=True)
class Field:
    name: str
    dtype: Any  # e.g. np.dtype("i4") or "U50"
    required: bool = True
    validator: Callable[[Any], None] | None = None  # raise on failure


@dataclass(frozen=True)
class Schema:
    fields: tuple[Field, ...]
    version: str
    dtype: np.dtype = field(init=False, repr=False)
    label_from: tuple[str, ...] | None = None
    label_required: bool = False
    label_builder: Callable[[tuple[str, ...]], str] | None = None
    label_normalizer: Callable[[str], str] | None = None
    _field_by_name: dict[str, Field] = field(init=False, repr=False)

    def __post_init__(self):
        object.__setattr__(self, "dtype", np.dtype(list(self.get_spec())))
        fb = {f.name: f for f in self.fields}

        if len(fb) != len(self.fields):
            raise ValueError("Duplicate field names in schema.")
        object.__setattr__(self, "_field_by_name", fb)

    def field(self, name: str) -> Field:
        try:
            return self._field_by_name[name]
        except KeyError:
            raise KeyError(f"Field {name!r} not in schema.") from None

    def required_names(self) -> tuple[str, ...]:
        return tuple(f.name for f in self.fields if f.required)

    def make_label(self, row: dict[str, Any]) -> str:
        """Build a label from `row` using `label_from`, `label_builder`, and `label_normalizer`.

        This is intended for Registry to call during add/update so callers don't need
        to pass the schema each time.
        """
        if self.label_from is None:
            raise ValueError("Schema.label_from is None; cannot build label")

        # Extract the requested fields and stringify.
        parts = tuple(str(row[k]) for k in self.label_from)

        # Default builder if none provided.
        builder = self.label_builder or (lambda xs: "|".join(xs))
        label = builder(parts)

        if self.label_normalizer is not None:
            label = self.label_normalizer(label)

        return label

    def get_spec(self):
        return tuple((field.name, field.dtype) for field in self.fields)

    def get_names(self):
        return tuple(field.name for field in self.fields)

    def schema_pairs(self):
        return [(k, str(v)) for k, v in self.get_spec()]

    def get_payload(self):
        return json.dumps(
            [(k, str(v)) for k, v in self.get_spec()], separators=(",", ":")
        )

    def json_pairs(self):
        return self.get_payload()

    def hash(self):
        payload = self.get_payload().encode("utf-8")
        return hashlib.sha256(payload).hexdigest()
