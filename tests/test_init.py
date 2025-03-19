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


def test_dnastream_str(dnastream_obj):
    """
    Test the __str__ method of DNAStream to ensure it runs correctly.
    """

    # Call __str__ method
    output = str(dnastream_obj)

    # Ensure the output is a non-empty string
    assert isinstance(output, str), "__str__ method did not return a string"
    assert len(output) > 0, "__str__ method returned an empty string"

    # Check for expected substrings in the output
    expected_substrings = [
        "DNAStream object",
        "SNVs",
        "samples",
        "HDF5 File:",
        "Patient:",
        "sex:",
    ]

    for substring in expected_substrings:
        assert substring in output, f"Missing expected text: {substring}"
