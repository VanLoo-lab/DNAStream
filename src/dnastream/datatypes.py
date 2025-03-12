# Description: This file contains the definition of the custom datatypes used in the DNAStream package.
import h5py
import numpy as np

STR_DTYPE = h5py.string_dtype(encoding="utf-8")
# datatype for metadata for SNV index
SNV_DTYPE = np.dtype(
    [
        ("label", STR_DTYPE),
        ("chrom", "S5"),  # Fixed-length string (5 characters max)
        ("pos", "i8"),
        ("end_pos", "i8"),  # Integer
        ("ref_allele", "S1"),  # Single-character string
        ("alt_allele", "S1"),  # Single-character string
        ("hugo", "S15"),  # Fixed-length string (15 characters max)
        ("gene", "S10"),  # Fixed-length string (10 characters ma)
    ]
)

# datatype to store metadata for sample index
SAMPLE_DTYPE = np.dtype(
    [
        ("label", STR_DTYPE),  # sample description
        ("patient", "S10"),  # Fixed-length string (10 characters max)
        ("source", "S10"),  # Fixed-length string (10 characters max)
        ("location", "S15"),  # Fixed-length string (15 characters max)
        ("file", STR_DTYPE),  #
    ]
)

# datatype to for the index modification logs
LOG_DTYPE = np.dtype(
    [
        ("timestamp", STR_DTYPE),
        ("number", "i8"),
        ("index_size_before", "i8"),
        ("index_size_after", "i8"),
        ("operation", "S15"),
        ("user", "S15"),
        ("hostname", "S15"),
        ("file", STR_DTYPE),
    ]
)

# datatype for tree edge list
EDGE_LIST_DTYPE = np.dtype([("parent", "i4"), ("child", "i4")])

# datatype for storing tree metadata
TREE_DTYPE = np.dtype(
    [
        ("label", "S15"),
        ("method", "S15"),
        ("score", "f8"),
        ("rank", "i4"),
        ("file", STR_DTYPE),
    ]
)

VLEN_EDGE_DTYPE = h5py.special_dtype(vlen=EDGE_LIST_DTYPE)

# datatype for the dataset track log
DATASET_LOG_DTYPE = np.dtype(
    [
        ("timestamp", STR_DTYPE),  # Timestamp
        ("dataset", STR_DTYPE),  # Dataset path modified
        ("operation", "S15"),  # Operation type: "add", "update", "delete"
        ("user", "S15"),  # User who performed the modification
        ("hostname", "S15"),  # Host machine
        ("file", STR_DTYPE),  # source file
    ]
)

SEGMENT_LABEL_DTYPE = np.dtype([("chrom", "S5"), ("start", "i8"), ("end", "i8")])

ALLELE_SPECIFIC_CN_DTYPE = np.dtype([("x", "i4"), ("y", "i4"), ("proportion", "f4")])
