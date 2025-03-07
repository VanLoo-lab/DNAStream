import pytest
import os
import h5py
from dnastream import DNAStream

# Test data paths (modify as needed)
RC_PATH = "/rsrch6/home/genetics/vanloolab/llweber/MPNST/scdna/read_counts"
MAF_FILES = [f"../data_summary/WGS_MUTATION/Somatic/2outof3_SNV/GEM2.2_PT_{i}_SNVs_2outof3.maf" for i in range(2, 6)]
READ_COUNT_FILE = f"{RC_PATH}/GEM2.2.csv"

@pytest.fixture
def temp_h5_file():
    """Creates a temporary HDF5 file for testing and cleans up afterward."""
    filename = "test_temp.h5"
    yield filename  # Test function runs here
    if os.path.exists(filename):
        os.remove(filename)

def test_dnastream_initialization(temp_h5_file):
    """Test DNAStream object initialization."""
    ds = DNAStream(filename=temp_h5_file, verbose=True)
    assert ds.file is not None
    assert isinstance(ds.file, h5py.File)
    ds.close()

def test_add_maf_files(temp_h5_file):
    """Test adding MAF files."""
    ds = DNAStream(filename=temp_h5_file, verbose=True)
    indices = ds.add_maf_files(MAF_FILES)
    assert indices is None

    ds.close()

def test_add_read_counts(temp_h5_file):
    """Test adding read counts."""
    ds = DNAStream(filename=temp_h5_file, verbose=True)
    ds.add_read_counts(READ_COUNT_FILE, source="scdna")
    snv_log = ds.get_snv_log()
    sample_log = ds.get_sample_log()
    
    assert not snv_log.empty
    assert not sample_log.empty
    ds.close()

def test_log_retrieval(temp_h5_file):
    """Test that logs are correctly retrieved after modifications."""
    ds = DNAStream(filename=temp_h5_file, verbose=True)
    ds.add_maf_files(MAF_FILES)
    ds.add_read_counts(READ_COUNT_FILE, source="scdna")
    
    snv_log = ds.get_snv_log()
    sample_log = ds.get_sample_log()
    
    assert snv_log is not None
    assert sample_log is not None
    assert len(snv_log) > 0
    assert len(sample_log) > 0
    ds.close()

def test_dnastream_cleanup(temp_h5_file):
    """Ensure the HDF5 file closes properly after operations."""
    ds = DNAStream(filename=temp_h5_file, verbose=True)
    ds.close()
    assert ds.file.id.valid == 0  # Check that the file is closed