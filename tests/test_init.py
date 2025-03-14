import h5py
from dnastream import DNAStream


def test_dnastream_initialization(temp_h5_file):
    """Test DNAStream object initialization."""
    dnastream_obj = DNAStream(filename=temp_h5_file, verbose=True)
    assert dnastream_obj.file is not None
    assert isinstance(dnastream_obj.file, h5py.File)
    dnastream_obj.close()


def test_dnastream_cleanup(dnastream_obj):
    """Ensure the HDF5 file closes properly after operations."""
    dnastream_obj.close()
    assert not dnastream_obj.file.id.valid  # File should be closed
