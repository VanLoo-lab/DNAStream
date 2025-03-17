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
    "metadata": {
        "dtype": SEGMENT_LABEL_DTYPE,
        "shape": (0,),
        "maxshape": (None,),
        "chunks": (100,),
    },
    "labels": LABEL_DICT,
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
    "index": {
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
            "tracked_tables": [("read_counts/variant", 0), ("read_counts/total", 0)],
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
            "tracked_tables": [
                (f"{ctab}/{ctype}", 1)
                for ctab in ["read_counts", "allele_counts"]
                for ctype in ["variant", "total"]
            ]
            + [
                (f"copy_numbers/{s}/{ctab}", 1)
                for ctab in ["profile", "logr", "baf"]
                for s in MODALITIES
            ],
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
            "tracked_tables": [
                ("allele_counts/variant", 0),
                ("allele_counts/total", 0),
            ],
        },
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
            "tracked_tables": [("tree/SNV_trees/trees", 0)],
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
            "tracked_tables": [("tree/CNA_trees/trees", 0)],
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
    "allele_counts": {
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

# add tracking for local segment index along axis 0
for s in MODALITIES:
    SCHEMA["copy_numbers"][s]["tracked_tables"] = []
    for tab in ["profile", "logr", "baf"]:
        SCHEMA["copy_numbers"][s]["tracked_tables"].append(
            (f"copy_numbers/{s}/{tab}", 0)
        )


def get_schema_value(path, key, default=None):
    """
    Retrieve a value from SCHEMA based on a given HDF5-like path.

    Parameters
    ----------
    path : str
        The hierarchical path in the SCHEMA dictionary (e.g., "index/SNV").
    key : str
        The specific key to retrieve within the resolved SCHEMA dictionary.
    default : any, optional
        The default value to return if the key is not found (default is None).

    Returns
    -------
    any
        The value associated with the key in SCHEMA, or the default if not found.
    """
    parts = path.split("/")
    schema_section = SCHEMA

    for part in parts:
        schema_section = schema_section.get(part, {})
        if not isinstance(schema_section, dict):
            return default  # Stop if we hit a non-dict value

    return schema_section.get(key, default)
