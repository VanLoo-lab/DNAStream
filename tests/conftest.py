import os
import pytest
import h5py
from dnastream import DNAStream
from dnastream.registry import Registry
from dnastream.schema import Schema, Field
import numpy as np


STR_DTYPE = h5py.string_dtype("utf-8")


REGISTRY_SPINE = (
    Field("id", STR_DTYPE, True, None),  # UUIDv4 string
    Field("label", STR_DTYPE, True, None),  # user-facing unique key
    Field("idx", np.int64, True, None),
    Field("active", np.bool_, True, None),
    Field("created_at", STR_DTYPE, True, None),  # ISO8601 Z
    Field("created_by", STR_DTYPE, True, None),
    Field("modified_at", STR_DTYPE, True, None),  # ISO8601 Z
    Field("modified_by", STR_DTYPE, True, None),
)


@pytest.fixture
def temp_h5_file(tmp_path):
    """Path to a temporary HDF5 file for testing."""
    return tmp_path / "test.h5"


@pytest.fixture
def temp_h5_handle(temp_h5_file):
    handle = h5py.File(temp_h5_file, "w")
    try:
        yield handle
    finally:
        handle.close()


@pytest.fixture
def temp_data_schema():
    temp_spec = (
        Field(name="variable", dtype=h5py.string_dtype("utf-8"), required=True),
        Field(name="value", dtype=np.int64, required=False),
    )
    return Schema(
        fields=temp_spec,
        version="1.0.1",
        label_from=("variable",),
        label_required=True,
        label_builder=lambda variable: (
            variable.lower() if isinstance(variable, str) else ""
        ),
        label_normalizer=lambda x: str(x),
    )


@pytest.fixture
def temp_registry_schema():
    temp_fields = REGISTRY_SPINE + (
        Field("variable", h5py.string_dtype("utf-8"), required=True),
        Field("value", np.int64),
    )

    return Schema(
        fields=temp_fields,
        version="1.0.1",
        label_from=("variable",),
        label_required=True,
        label_builder=lambda x: str(x[0]).lower(),
        label_normalizer=lambda x: str(x),
    )


@pytest.fixture
def temp_data_rows():
    temp_data = [{"variable": f"var{i}", "value": i} for i in range(5)]
    return temp_data


@pytest.fixture
def temp_registry(registry_obj, temp_data_schema):

    return registry_obj.create(temp_data_schema)


@pytest.fixture
def dnastream_obj(temp_h5_file):
    """Provides a connected DNAStream instance for testing."""
    ds = DNAStream(str(temp_h5_file), mode="x")
    ds.connect()
    try:
        yield ds
    finally:
        ds.close()


@pytest.fixture
def registry_obj(temp_h5_handle, temp_registry_schema):
    grp = temp_h5_handle.require_group("registry")
    reg = Registry(grp, "test")
    reg.create(temp_registry_schema)
    return reg


# @pytest.fixture
# def temp_h5_stream(temp_h5_file):
#     """Fixture to create a temporary HDF5 file for testing."""

#     with h5py.File(temp_h5_file, "a") as f:
#         f.create_dataset(
#             "read_counts",
#             shape=(0, 0),
#             dtype="i",
#             maxshape=(None, None),
#             chunks=(1, 5000),
#         )
#         yield f
#         f.close()


# @pytest.fixture
# def base_index(temp_h5_stream):
#     """Fixture to create a BaseIndex instance."""

#     yield LocalIndex(
#         temp_h5_stream,
#         name="index/SNV",
#         metadata_dtype=VariantMetadata.get_dtype(),
#         verbose=True,
#         tracked_tables=[("read_counts", 0)],
#     )


# @pytest.fixture
# def tree_file():
#     """Provides the path to the tree file for testing."""
#     return os.path.join(TEST_DATA_DIR, "trees.txt")


# @pytest.fixture
# def battenberg_file():
#     """Provides the path to the tree file for testing."""
#     return os.path.join(TEST_DATA_DIR, "battenberg.tsv")


# @pytest.fixture
# def global_index(temp_h5_stream):
#     """Fixture to create a BaseIndex instance."""

#     yield GlobalIndex(
#         temp_h5_stream,
#         name="index/SNV",
#         metadata_dtype=VariantMetadata.get_dtype(),
#         tracked_tables=[("read_counts", 0)],
#         verbose=True,
#     )


# @pytest.fixture
# def pyclone_file():
#     """Fixture to provide pyclone file for testing."""
#     return os.path.join(TEST_DATA_DIR, "pyclone.tsv")


# @pytest.fixture
# def ascat_total_file():
#     """Fixture to provide pyclone file for testing."""
#     return os.path.join(TEST_DATA_DIR, "ASCAT.scprofile.txt")

# @pytest.fixture
# def ascat_as_file():
#     """Fixture to provide ASCATsc allele-specific file for testing."""
#     return os.path.join(TEST_DATA_DIR, "AS_ascat_profile.txt")


# @pytest.fixture
# def read_count_file():
#     """Fixture to provide pyclone file for testing."""
#     return os.path.join(TEST_DATA_DIR, "read_counts.csv")


# @pytest.fixture
# def sample_metadata_file():
#     """Fixture to provide pyclone file for testing."""
#     return os.path.join(TEST_DATA_DIR, "sample_metadata.csv")


# @pytest.fixture
# def maf_file():
#     """Provides the path to the maf file for testing."""
#     return os.path.join(TEST_DATA_DIR, "test.maf")
