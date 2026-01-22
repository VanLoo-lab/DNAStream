import pytest
from dnastream.schema import Schema, Field
import numpy as np


def test_raises_on_duplicate_field_names():
    """
    Test that duplicate field names raises an error on initialization.
    """

    with pytest.raises(ValueError, match="Duplicate field names"):
        Schema(fields=(Field("x", "S10"), Field("x", "S10")), version="1.0.0")


def test_raises_on_invalid_label_from():
    """
    Test that invalid combinations of label_required and label_from raise an error on init
    and that label_from only contains fields in the schema.
    """
    fields = (Field("x", "S10"), Field("y", np.int64))
    with pytest.raises(ValueError):
        Schema(fields=fields, version="1.0.0", label_required=True, label_from=None)
        Schema(fields=fields, version="1.0.0", label_required=True, label_from=("z",))


def test_dtype_field_constructed_properly():
    """
    Test that that an np.dtype is correctly constructed post initialization.
    """
    fields = (Field("x", "S10"), Field("y", np.int64))
    my_schema = Schema(fields=fields, version="1.0.0")

    for name, fld in zip(my_schema.dtype.names, fields):
        assert name == fld.name
        assert my_schema.dtype[fld.name] == np.dtype(fld.dtype)


def test_fields_function():
    """
    Test the `field` function correctly returns the Field from its name and raises KeyError when
    field is not present in the schema.
    """
    test_field = Field("x", "S10")
    my_schema = Schema(fields=(test_field, Field("y", "S10")), version="1.0.0")

    with pytest.raises(KeyError):
        my_schema.field("foo")

    my_field = my_schema.field("x")
    assert test_field == my_field


def test_hash_function_comparison():
    """
    test that two schemas that contain the same set of fields are identical.
    The order of the fields doesn't matter.
    """
    fields = (Field("x", "S10"), Field("y", np.int64))
    my_schema = Schema(fields=fields, version="1.0.0")
    my_schema2 = Schema(fields=fields, version="1.0.0")
    assert my_schema.version == my_schema2.version
    assert my_schema == my_schema2
    assert my_schema.hash() == my_schema2.hash()
    assert my_schema.json_pairs() == my_schema2.json_pairs()

    # field order
    my_schema3 = Schema(
        fields=(
            Field("y", np.int64),
            Field("x", "S10"),
        ),
        version="1.0.0",
    )
    assert my_schema3 != my_schema
    assert my_schema3.hash() == my_schema.hash()
    assert my_schema.json_pairs() == my_schema3.json_pairs()


def test_make_label():
    """
    Test make label generates correct label and raises ValueError when label_from isn't provided.
    """
    fields = (Field("x", "S10"), Field("y", np.int64))
    row = {"x": "X", "y": "Y"}
    my_schema = Schema(fields=fields, version="1.0.0")

    with pytest.raises(ValueError):
        lab = my_schema.make_label(row)
    # assert my_schema.make_label(row) is None

    my_schema = Schema(
        fields=fields,
        version="1.0.0",
        label_from=("x", "y"),
        label_required=True,
        label_builder=lambda x: "|".join(x),
    )
    row = {"x": "X", "y": "Y"}
    assert "X|Y" == my_schema.make_label(row)

    my_schema = Schema(
        fields=fields,
        version="1.0.0",
        label_from=("x", "y"),
        label_required=True,
        label_builder=lambda x: "|".join(x),
        label_normalizer=lambda x: x.lower(),
    )
    assert "x|y" == my_schema.make_label(row)


def test_required_names():
    """
    Test to ensure required names returns the names of all required fields
    """
    fields = (
        Field("x", "S10", required=True),
        Field("y", np.int64, required=True),
        Field("z", np.int16, required=False),
    )

    my_schema = Schema(fields=fields, version="1.0.0")
    required = {"x", "y"}
    assert set(my_schema.required_names()) == required
