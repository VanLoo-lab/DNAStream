import pytest
from dnastream.schema import Schema, Field
import numpy as np


def test_raises_on_duplicate_field_names():

    with pytest.raises(ValueError, match="Duplicate field names"):
        Schema(fields=(Field("x", "S10"), Field("x", "S10")), version="1.0.0")


def test_dtype_field_constructed_properly():
    fields = (Field("x", "S10"), Field("y", np.int64))
    my_schema = Schema(fields=fields, version="1.0.0")

    for name, fld in zip(my_schema.dtype.names, fields):
        assert name == fld.name
        assert my_schema.dtype[fld.name] == np.dtype(fld.dtype)


def test_fields_function():
    test_field = Field("x", "S10")
    my_schema = Schema(fields=(test_field, Field("y", "S10")), version="1.0.0")

    with pytest.raises(KeyError):
        my_schema.field("foo")

    my_field = my_schema.field("x")
    assert test_field == my_field
