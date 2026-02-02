import pytest
import h5py
from tests.helpers.h5dataset import _TestH5Dataset
from dataclasses import replace
from dnastream.schema import Field


"""
Tests the internal low-level H5Dataset class from which
more specific and specialized classes inherit from, such as 
Registry and Measurement, Provenance 
"""


def test_h5dataset_init(temp_h5_handle):
    """
    Given a handle to an h5 file, a group and a dataset name,
    the H5Dataset object should be initialized with proper path set
    but the dataset not yet created in the h5 file.
    """
    parent = temp_h5_handle.require_group("data")
    dt = _TestH5Dataset(parent, "my_data")

    assert dt.path == "/data/my_data"
    assert not dt.exists()


def test_h5dataset_create(temp_h5_handle, temp_data_schema):
    """
    Given a handle to an h5 file and schema, the dataset should be
    created with the schema specific attributes and accessible. `schema_version` and `schema_hash`
    specificed as dataset attributes.  The dtype of the dataset must also match the schema.
    """
    parent = temp_h5_handle.require_group("data")
    dt = _TestH5Dataset(parent, "my_data")

    ds = dt.create(temp_data_schema)
    assert isinstance(ds, h5py.Dataset)
    handle = dt.open(temp_data_schema, strict=True)
    assert handle
    assert handle.name == "/data/my_data"
    assert handle.dtype.descr == temp_data_schema.dtype.descr
    assert handle.attrs["schema_version"] == str(temp_data_schema.version)
    assert handle.attrs["schema_hash"] == str(temp_data_schema.hash())
    assert "schema_pairs_json" in handle.attrs


def test_h5dataset_open_raises_on_schema_version_mismatch(
    temp_h5_handle, temp_data_schema
):
    """
    Given a handle to an h5 file and schema version mismatches are caught and
    a ValueError is raised when strict is True.
    """
    parent = temp_h5_handle.require_group("data")
    dt = _TestH5Dataset(parent, "my_data")
    dt.create(temp_data_schema)

    bad = replace(temp_data_schema, version="9.9.9")

    with pytest.raises(ValueError, match="schema_version mismatch"):
        dt.open(bad, strict=True)


def test_h5dataset_open_raises_on_schema_hash_mismatch(
    temp_h5_handle, temp_data_schema
):
    """
    Given a handle to an h5 file and schema hash mismatches are caught and
    a ValueError is raised when strict is true.
    """
    parent = temp_h5_handle.require_group("data")
    dt = _TestH5Dataset(parent, "my_data")
    dt.create(temp_data_schema)

    bad = replace(
        temp_data_schema,
        fields=(
            Field(name="foo", dtype="S10"),
            Field("variable", h5py.string_dtype("utf-8"), required=True),
        ),
    )

    with pytest.raises(ValueError, match="schema_hash mismatch"):
        dt.open(bad, strict=True)


def test_h5dataset_open_warns_when_not_strict(temp_h5_handle, temp_data_schema):
    """
    Given an existing dataset with a stored schema identity,
    when an expected schema mismatches and strict=False,
    then open() should warn and still return the dataset handle.
    """
    parent = temp_h5_handle.require_group("data")
    dt = _TestH5Dataset(parent, "my_data")
    dt.create(temp_data_schema)

    bad = replace(
        temp_data_schema,
        fields=(
            Field(name="foo", dtype="S10"),
            Field("variable", h5py.string_dtype("utf-8"), required=True),
        ),
    )

    with pytest.warns(UserWarning):
        dt.open(bad, strict=False)


def test_h5dataset_open_missing_raises(temp_h5_handle, temp_data_schema):
    """
    Given an existing dataset with a stored schema identity,
    when the dataset is attempted to be opened but was not yet created,
    then open() should raise a RuntimeError.
    """
    parent = temp_h5_handle.require_group("data")
    dt = _TestH5Dataset(parent, "missing")
    with pytest.raises(RuntimeError, match="does not exist"):
        dt.open(temp_data_schema)


def test_h5dataset_create_refuses_if_exists(temp_h5_handle, temp_data_schema):
    """
    Given an existing dataset with a stored schema identity,
    when the dataset is attempted to be createed multiple times
    then create() throws a RunTimeError
    """
    parent = temp_h5_handle.require_group("data")
    dt = _TestH5Dataset(parent, "my_data")
    dt.create(temp_data_schema)

    with pytest.raises(RuntimeError, match="already exists"):
        dt.create(temp_data_schema)


def test_h5dataset_create_refuses_dtype(temp_h5_handle, temp_data_schema):
    """
    Given an existing dataset with a stored schema identity,
    when the kwargs contain dtype then create() throws a TypeError
    """
    parent = temp_h5_handle.require_group("data")
    dt = _TestH5Dataset(parent, "my_data")
    with pytest.raises(TypeError):
        dt.create(temp_data_schema, dtype="i4")
