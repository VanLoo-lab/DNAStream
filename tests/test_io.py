def test_add_pyclone(dnastream_obj, pyclone_file):
    dnastream_obj.parse_pyclone_file(pyclone_file)
    assert (
        dnastream_obj.get_snv_size() > 0
    ), "No SNVs were added to the DNAStream object"
    assert (
        dnastream_obj.get_sample_size() > 0
    ), "No samples were added to the DNAStream object"
    cluster_dict = dnastream_obj.get_snv_clusters()
    assert len(cluster_dict) > 0, "No clusters were added to the DNAStream object"
    assert cluster_dict["m0"] == 0, "SNV m0 should be in cluster 0"
    assert cluster_dict["m29"] == 3, "SNV m1 should be in cluster 1"


def test_add_battenberg(dnastream_obj, battenberg_file):
    dnastream_obj.parse_battenberg_file(battenberg_file, "my_test_sample")
    assert (
        dnastream_obj.get_snv_size() == 0
    ), "No SNVs were added to the DNAStream object"
    assert (
        dnastream_obj.get_sample_size() == 1
    ), "One sample should be added to the DNAStream object"
    assert (
        dnastream_obj.get_copy_numbers_bulk_size() > 0
    ), "No segments were added to the copy_numbers_bulk local index!"

    indices = dnastream_obj.indices_by_sample_label(["my_test_sample"])
    assert len(indices) == 1, "One sample should be added to the DNAStream object"
    assert indices[0] == 0, "Sample index should be 0"
    sample_idx = indices[0]
    seg_indices = dnastream_obj.indices_by_copy_numbers_bulk_label(["1:100000:200000"])
    assert len(seg_indices) == 1, "One segment should be added to the DNAStream object"
    assert seg_indices[0] == 0, "Segment index should be 0"
    assert (
        dnastream_obj["copy_numbers/bulk/logr"][seg_indices[0], sample_idx] == 0.2
    ), "logr for 0,0 should be 0.2"
    assert (
        dnastream_obj["copy_numbers/bulk/baf"][seg_indices[0], sample_idx] == 0.45
    ), "baf for 0,0 should be 0.2"


def test_add_ascatsc(dnastream_obj, ascat_total_file):
    dnastream_obj.parse_ascat_sc_total_copy_number_file(
        ascat_total_file, sample_label="my_test_sample"
    )
    assert (
        dnastream_obj.get_snv_size() == 0
    ), "No SNVs were added to the DNAStream object"
    assert (
        dnastream_obj.get_sample_size() == 1
    ), "One sample should be added to the DNAStream object"
    assert (
        dnastream_obj.get_copy_numbers_scdna_size() > 0
    ), "No segments were added to the copy_numbers_scdna local index!"

    indices = dnastream_obj.indices_by_sample_label(["my_test_sample"])
    assert len(indices) == 1, "One sample should be added to the DNAStream object"
    assert indices[0] == 0, "Sample index should be 0"
    sample_idx = indices[0]
    seg_indices = dnastream_obj.indices_by_copy_numbers_scdna_label(
        ["chr1:617509:8524993"]
    )
    assert len(seg_indices) == 1, "One segment should be added to the DNAStream object"
    assert seg_indices[0] == 0, "Segment index should be 0"
    assert (
        dnastream_obj["copy_numbers/scdna/logr"][seg_indices[0], sample_idx]
        == -0.716912189895419
    ), "logr for 0,0 should be 0.2"


def test_add_read_counts(dnastream_obj, read_count_file):
    """Test adding read counts."""

    dnastream_obj.add_read_counts(
        read_count_file, source="scdna", columns={"cell": "sample"}
    )
    snv_log = dnastream_obj.get_snv_log()
    sample_log = dnastream_obj.get_sample_log()

    assert not snv_log.empty
    assert not sample_log.empty
