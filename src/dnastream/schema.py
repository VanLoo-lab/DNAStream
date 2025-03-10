import h5py
from .datatypes import (
    LOG_DTYPE,
    SNV_DTYPE,
    SAMPLE_DTYPE,
    VLEN_EDGE_DTYPE,
    DATASET_LOG_DTYPE,
    TREE_DTYPE,
)


MODALITIES = ["bulk", "lcm", "scdna"]

STRUCT_ARRAYS = ["log", "data"]

META_TABLES = ["data", "label", "cluster", "index", "log"]


READ_COUNTS = {
    f"read_counts/{m}/{c}": "i" for c in ["variant", "total"] for m in MODALITIES
}

SCHEMA = {
    "metadata": {"log": DATASET_LOG_DTYPE},
    "SNV": {
        "label": h5py.string_dtype(encoding="utf-8"),
        "cluster": "i8",
        "data": SNV_DTYPE,
        "index": h5py.string_dtype("utf-8"),
        "log": LOG_DTYPE,
    },
    "sample": {
        "label": h5py.string_dtype(encoding="utf-8"),
        "cluster": "i8",
        "data": SAMPLE_DTYPE,
        "index": h5py.string_dtype("utf-8"),
        "log": LOG_DTYPE,
    },
    "tree": {
        "SNV_trees": {
            "trees": VLEN_EDGE_DTYPE,
            "data": TREE_DTYPE,
            "index": h5py.string_dtype("utf-8"),
        },
        "CNA_trees": {
            "trees": VLEN_EDGE_DTYPE,
            "data": TREE_DTYPE,
            "index": h5py.string_dtype("utf-8"),
        },
        "clonal_trees": {
            "trees": VLEN_EDGE_DTYPE,
            "data": TREE_DTYPE,
            "index": h5py.string_dtype("utf-8"),
        },
    },
}
