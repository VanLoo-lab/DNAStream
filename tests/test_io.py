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
