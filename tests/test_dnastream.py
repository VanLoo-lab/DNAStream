import pytest
import os
import h5py
from dnastream import DNAStream

# Test data paths (modify as needed)
RC_PATH = "/rsrch6/home/genetics/vanloolab/llweber/MPNST/scdna/read_counts"
MAF_FILES = [
    f"/rsrch6/home/genetics/vanloolab/llweber/data_summary/WGS_MUTATION/Somatic/2outof3_SNV/GEM2.2_PT_{i}_SNVs_2outof3.maf"
    for i in range(2, 6)
]
READ_COUNT_FILE = f"{RC_PATH}/GEM2.2.csv"
CONIPHER = "/rsrch6/home/genetics/vanloolab/llweber/MPNST/tree_building/trees/GEM2.2/conipher/trees.txt"
SAPLING = "/rsrch6/home/genetics/vanloolab/llweber/MPNST/tree_building/trees/GEM2.2/sapling/trees.txt"


# def test_add_maf_files(temp_h5_file):
#     """Test adding MAF files."""
#     ds = DNAStream(filename=temp_h5_file, verbose=True)
#     indices = ds.add_maf_files(MAF_FILES)
#     assert indices is None

#     ds.close()


# def test_add_read_counts(temp_h5_file):
#     """Test adding read counts."""
#     ds = DNAStream(filename=temp_h5_file, verbose=True)
#     ds.add_read_counts(READ_COUNT_FILE, source="scdna")
#     snv_log = ds.get_snv_log()
#     sample_log = ds.get_sample_log()

#     assert not snv_log.empty
#     assert not sample_log.empty
#     ds.close()


# def test_log_retrieval(temp_h5_file):
#     """Test that logs are correctly retrieved after modifications."""
#     ds = DNAStream(filename=temp_h5_file, verbose=True)
#     ds.add_maf_files(MAF_FILES)
#     ds.add_read_counts(READ_COUNT_FILE, source="scdna")
#     snv_log = ds.get_snv_log()
#     sample_log = ds.get_sample_log()

#     assert snv_log is not None
#     assert sample_log is not None
#     assert len(snv_log) > 0
#     assert len(sample_log) > 0
#     ds.close()
