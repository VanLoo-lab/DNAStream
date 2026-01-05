import os
import pytest
import h5py
from dnastream import DNAStream
from dnastream.registry import Registry
import numpy as np
import hashlib
import json


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
    temp_spec = (("variable", h5py.string_dtype("utf-8")), ("value", np.int64))

    payload = json.dumps(
        [(k, str(v)) for k, v in temp_spec],
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")

    temp_schema = {
        "spec": temp_spec,
        "dtype": np.dtype(list(temp_spec)),
        "columns": [k for k, _ in temp_spec],
        "schema_pairs": [(k, str(v)) for k, v in temp_spec],
        "schema_hash": hashlib.sha256(payload).hexdigest(),
        "schema_version": "1.0.1",
        "label_from": "variable",
        "label_normalizer": lambda x: x.lower() if isinstance(x, str) else "",
        "label_required": False,
    }
    return temp_schema


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
def registry_obj(temp_h5_handle):
    reg = Registry(temp_h5_handle, "")


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
