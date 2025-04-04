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
