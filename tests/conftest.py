import os
import pytest
import h5py
from dnastream import DNAStream
from dnastream.index_manager import LocalIndex, GlobalIndex
from dnastream.datatypes import SNV_DTYPE, STR_DTYPE

# Define test data directory
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture
def temp_h5_file():
    """Creates a temporary HDF5 file for testing and cleans up afterward."""
    filename = os.path.join(TEST_DATA_DIR, "test_temp.h5")
    yield filename  # Test function runs here
    if os.path.exists(filename):
        os.remove(filename)


@pytest.fixture
def temp_h5_stream(temp_h5_file):
    """Fixture to create a temporary HDF5 file for testing."""

    with h5py.File(temp_h5_file, "a") as f:
        f.create_dataset(
            "read_counts",
            shape=(0, 0),
            dtype="i",
            maxshape=(None, None),
            chunks=(1, 5000),
        )
        yield f
        f.close()


@pytest.fixture
def base_index(temp_h5_stream):
    """Fixture to create a BaseIndex instance."""

    yield LocalIndex(
        temp_h5_stream,
        name="index/SNV",
        metadata_dtype=SNV_DTYPE,
        verbose=True,
        tracked_tables=[("read_counts", 0)],
    )


@pytest.fixture
def dnastream_obj(temp_h5_file):
    """Provides a DNAStream instance for testing."""
    return DNAStream(temp_h5_file, initialize=True, verbose=False)


@pytest.fixture
def tree_file():
    """Provides the path to the tree file for testing."""
    return os.path.join(TEST_DATA_DIR, "trees.txt")


@pytest.fixture
def global_index(temp_h5_stream):
    """Fixture to create a BaseIndex instance."""

    yield GlobalIndex(
        temp_h5_stream,
        name="index/SNV",
        metadata_dtype=SNV_DTYPE,
        tracked_tables=[("read_counts", 0)],
        verbose=True,
    )


@pytest.fixture
def pyclone_file():
    """Fixture to provide pyclone file for testing."""
    return os.path.join(TEST_DATA_DIR, "pyclone.tsv")
