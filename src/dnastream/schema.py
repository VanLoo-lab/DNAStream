# Description: This file contains the schema for the DNAStream metadata model.
import h5py
from .datatypes import (
    STR_DTYPE,
    LOG_DTYPE,
    SNV_DTYPE,
    SAMPLE_DTYPE,
    VLEN_EDGE_DTYPE,
    DATASET_LOG_DTYPE,
    TREE_DTYPE,
    SEGMENT_LABEL_DTYPE,
    ALLELE_SPECIFIC_CN_DTYPE,
)


MODALITIES = ["bulk", "lcm", "scdna"]

STRUCT_ARRAYS = ["log", "metadata"]

META_TABLES = ["metadata", "cluster"]



INDEX_DICT = {
    "dtype": STR_DTYPE,
    "shape": (0,),
    "maxshape": (None,),
    "chunks": (100,),  # Chunking for efficient index expansion
}

LOG_DICT = {
    "dtype": LOG_DTYPE,
    "shape": (0,),
    "maxshape": (None,),
    "chunks": (100,),  # Log entries will be added incrementally
}

LABEL_DICT = {
    "dtype": STR_DTYPE,
    "shape": (0,),
    "maxshape": (None,),
    "chunks": (100,),
}

COPY_NUMBER_LAYER_DICT = {
    "profile": {
        "dtype": h5py.vlen_dtype(ALLELE_SPECIFIC_CN_DTYPE),
        "shape": (0, 0),
        "maxshape": (None, None),
        "chunks": (1, 100),  # Optimized for growing rows dynamically
    },
    "logr": {
        "dtype": "f8",
        "shape": (0, 0),
        "maxshape": (None, None),
        "chunks": (100, 100),
    },
    "baf": {
        "dtype": "f8",
        "shape": (0, 0),
        "maxshape": (None, None),
        "chunks": (100, 100),
    },
    "metadata":  {
            "dtype": SEGMENT_LABEL_DTYPE,
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
    "labels": LABEL_DICT
    # TODO: add metadata
}

SCHEMA = {
    "metadata": {
        "log": {
            "dtype": DATASET_LOG_DTYPE,
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        }
    },
    "index" : 
    {
        "SNV": {
        "cluster": {
            "dtype": "i8",
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
        "metadata": {
            "dtype": SNV_DTYPE,
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
        "labels": LABEL_DICT,
        "log": LOG_DICT,
    },
    "sample": {
        "label": {
            "dtype": STR_DTYPE,
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
        "cluster": {
            "dtype": "i8",
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
        "metadata": {
            "dtype": SAMPLE_DTYPE,
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
        "labels": LABEL_DICT,
        "log": LOG_DICT,
    },
    "SNP": {
        "label": {
            "dtype": STR_DTYPE,
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
        "cluster": {
            "dtype": "i8",
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
        "metadata": {
            "dtype": SNV_DTYPE,
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (100,),
        },
        "labels": LABEL_DICT,
        "log": LOG_DICT,
    }
    },
    "tree": {
        "SNV_trees": {
            "trees": {
                "dtype": VLEN_EDGE_DTYPE,
                "shape": (0,),
                "maxshape": (None,),
                "chunks": (50,),
            },
            "metadata": {
                "dtype": TREE_DTYPE,
                "shape": (0,),
                "maxshape": (None,),
                "chunks": (50,),
            },
            "labels": LABEL_DICT,
        },
        "CNA_trees": {
            "trees": {
                "dtype": STR_DTYPE,
                "shape": (0,),
                "maxshape": (None,),
                "chunks": (50,),
            },
            "metadata": {
                "dtype": TREE_DTYPE,
                "shape": (0,),
                "maxshape": (None,),
                "chunks": (50,),
            },
            "labels": LABEL_DICT,
        },
    },
    "copy_numbers": {s: COPY_NUMBER_LAYER_DICT for s in MODALITIES},
    "read_counts": {
        "variant": {
            "dtype": "i",
            "shape": (0, 0),
            "maxshape": (None, None),
            "chunks": (1, 5000),  # Optimized for row-wise updates
        },
        "total": {
            "dtype": "i",
            "shape": (0, 0),
            "maxshape": (None, None),
            "chunks": (1, 5000),
        },
        "log": LOG_DICT,
    },
}
