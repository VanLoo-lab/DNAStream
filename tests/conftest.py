import os
import pytest
from dnastream import DNAStream

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
def dnastream_obj(temp_h5_file):
    """Provides a DNAStream instance for testing."""
    return DNAStream(temp_h5_file, initialize=True, verbose=False)


@pytest.fixture
def tree_file():
    """Provides the path to the tree file for testing."""
    return os.path.join(TEST_DATA_DIR, "trees.txt")
