import h5py
import numpy as np


def test_base_index_init(base_index):
    """Test if the BaseIndex initializes correctly with empty datasets."""
    assert base_index.size() == 0, "Index should be empty!"

    expected_datasets = {"index/SNV/metadata", "index/SNV/labels"}
    for dat in expected_datasets:
        assert dat in base_index.file, f"{dat} not found in HDF5 file!"


def test_add_to_base_index(base_index):
    """Test adding entries to the BaseIndex."""
    base_index.add(labels=["SNV1", "SNV2"])

    assert base_index.size() == 2, "Index should have 2 elements!"
    assert (
        len(base_index._metadata_cache) == 2
    ), "Metadata cache should have 2 elements!"


def test_global_index_init(global_index):
    """Test if the GlobalIndex initializes correctly with empty datasets."""
    assert global_index.size() == 0, "Index should be empty!"

    expected_datasets = {
        "index/SNV/metadata",
        "index/SNV/labels",
        "index/SNV/log",
        "index/SNV/cluster",
    }
    for dat in expected_datasets:
        assert dat in global_index.file, f"{dat} not found in HDF5 file!"

    assert global_index.get_log().shape[0] == 1, "Log should be empty!"


def test_add_to_global_index(global_index):
    """Test adding entries to the BaseIndex."""
    global_index.print_log()
    global_index.add(labels=["SNV1", "SNV2"])
    global_index.print_log()

    assert global_index.size() == 2, "Index should have 2 elements!"
    # assert len(global_index._metadata_cache) == 2, "Metadata cache should have 2 elements!"
    assert global_index.get_log().shape[0] == 2, "Log should have 1 entry!"


def test_save_index(global_index):
    """Test saving the index to disk."""
    global_index.add(labels=["SNV1", "SNV2"])
    global_index.save_index()
    assert (
        global_index.last_saved_timestamp is not None
    ), "Timestamp should be set after saving!"
