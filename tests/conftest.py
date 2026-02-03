import os
import pytest
import h5py
from dnastream import DNAStream
from dnastream.registry import Registry
from dnastream.schema import Schema, Field
import numpy as np
from textwrap import dedent


STR_DTYPE = h5py.string_dtype("utf-8")


REGISTRY_SPINE = (
    Field("id", STR_DTYPE, True, None),  # UUIDv4 string
    Field("label", STR_DTYPE, True, None),  # user-facing unique key
    Field("idx", np.int64, True, None),
    Field("active", np.bool_, True, None),
    Field("created_at", STR_DTYPE, True, None),  # ISO8601 Z
    Field("created_by", STR_DTYPE, True, None),
    Field("modified_at", STR_DTYPE, True, None),  # ISO8601 Z
    Field("modified_by", STR_DTYPE, True, None),
)
_SPINE = {
    "id",
    "label",
    "idx",
    "active",
    "created_at",
    "created_by",
    "modified_at",
    "modified_by",
}


@pytest.fixture
def temp_h5_file(tmp_path):
    """Path to a temporary HDF5 file for testing."""
    return tmp_path / "test.h5"


@pytest.fixture
def temp_csv_file(tmp_path):
    """Path to a temporary HDF5 file for testing."""
    return tmp_path / "test.csv"


@pytest.fixture
def temp_h5_handle(temp_h5_file):
    handle = h5py.File(temp_h5_file, "w")
    try:
        yield handle
    finally:
        handle.close()


@pytest.fixture
def temp_data_schema():
    temp_spec = (
        Field(name="variable", dtype=h5py.string_dtype("utf-8"), required=True),
        Field(name="value", dtype=np.int64, required=False),
    )
    return Schema(
        fields=temp_spec,
        version="1.0.1",
        label_from=("variable",),
        label_required=True,
        label_builder=lambda variable: (
            variable.lower() if isinstance(variable, str) else ""
        ),
        label_normalizer=lambda x: str(x),
    )


@pytest.fixture
def temp_registry_schema():
    temp_fields = REGISTRY_SPINE + (
        Field("variable", h5py.string_dtype("utf-8"), required=True),
        Field("value", np.int64),
    )

    return Schema(
        fields=temp_fields,
        version="1.0.1",
        label_from=("variable",),
        label_required=True,
        label_builder=lambda x: str(x).lower(),
        label_normalizer=lambda x: str(x),
    )


@pytest.fixture
def temp_data_rows():
    temp_data = [{"variable": f"var{i}", "value": i} for i in range(5)]
    return temp_data


@pytest.fixture
def temp_registry(registry_obj, temp_data_schema):

    return registry_obj.create(temp_data_schema)


@pytest.fixture
def dnastream_obj(temp_h5_file):
    """Provides a connected DNAStream instance for testing."""
    ds = DNAStream.create(str(temp_h5_file))

    try:
        yield ds
    finally:
        ds.close()


@pytest.fixture
def registry_obj(temp_h5_handle, temp_registry_schema):
    grp = temp_h5_handle.require_group("registry")
    reg = Registry(grp, "test")
    reg.create(temp_registry_schema)
    return reg


@pytest.fixture
def temp_maf(tmp_path, name="test.maf", rows=None):
    path = tmp_path / name
    header = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Hugo_Symbol",
        "Entrez_Gene_Id",
        "Filter",
        "Variant_Classification",
        "Variant_Type",
        "dbSNP_RS",
    ]
    if rows is None:
        rows = [
            [
                "chr1",
                "123",
                "123",
                "A",
                "G",
                "TP53",
                "7157",
                "PASS",
                "Missense_Mutation",
                "SNP",
                "rs1",
            ],
            [
                "chr2",
                "456",
                "456",
                "C",
                "T",
                "EGFR",
                "1956",
                "PASS",
                "Missense_Mutation",
                "SNP",
                "rs2",
            ],
        ]

    lines = ["\t".join(header)] + ["\t".join(map(str, r)) for r in rows]
    path.write_text("\n".join(lines) + "\n")

    return path


def _dtype_kind(dt):
    # h5py.string_dtype ends up kind 'O' in numpy
    try:
        return np.dtype(dt).kind
    except Exception:
        return "O"


@pytest.fixture
def temp_sample_csv(dnastream_obj, tmp_path, name="samples.csv", rows=None):
    path = tmp_path / name

    # get the Schema object from the sample registry
    schema = getattr(dnastream_obj.sample, "schema", None) or getattr(
        dnastream_obj.sample, "_schema"
    )
    if schema is None:
        raise RuntimeError(
            "Couldn't find sample schema on dnastream_obj.sample (expected .schema or ._schema)"
        )

    # columns = all schema fields except registry spine
    header = [f.name for f in schema.fields if f.name not in _SPINE]

    # try to ensure the label_from fields exist and are set to something unique
    label_from = tuple(schema.label_from or ())

    if rows is None:
        rows = []
        for i in range(2):
            row = {}
            for f in schema.fields:
                if f.name in _SPINE:
                    continue
                kind = _dtype_kind(f.dtype)
                if kind in ("i", "u"):
                    row[f.name] = i
                elif kind == "f":
                    row[f.name] = float(i)
                elif kind == "b":
                    row[f.name] = True  # or set explicitly per-test
                else:
                    row[f.name] = f"{f.name}_{i}"

            # make label_from components nice + unique
            for j, key in enumerate(label_from):
                if key in row:
                    row[key] = f"{key}{i}"

            rows.append(row)

    # write CSV
    lines = [",".join(header)]
    for r in rows:
        lines.append(",".join(str(r.get(col, "")) for col in header))
    path.write_text("\n".join(lines) + "\n")

    return path


# @pytest.fixture
# def temp_h5_stream(temp_h5_file):
#     """Fixture to create a temporary HDF5 file for testing."""

#     with h5py.File(temp_h5_file, "a") as f:
#         f.create_dataset(
#             "read_counts",
#             shape=(0, 0),
#             dtype="i",
#             maxshape=(None, None),
#             chunks=(1, 5000),
#         )
#         yield f
#         f.close()


# @pytest.fixture
# def tree_file():
#     """Provides the path to the tree file for testing."""
#     return os.path.join(TEST_DATA_DIR, "trees.txt")


# @pytest.fixture
# def battenberg_file():
#     """Provides the path to the tree file for testing."""
#     return os.path.join(TEST_DATA_DIR, "battenberg.tsv")


# @pytest.fixture
# def global_index(temp_h5_stream):
#     """Fixture to create a BaseIndex instance."""

#     yield GlobalIndex(
#         temp_h5_stream,
#         name="index/SNV",
#         metadata_dtype=VariantMetadata.get_dtype(),
#         tracked_tables=[("read_counts", 0)],
#         verbose=True,
#     )


# @pytest.fixture
# def pyclone_file():
#     """Fixture to provide pyclone file for testing."""
#     return os.path.join(TEST_DATA_DIR, "pyclone.tsv")


# @pytest.fixture
# def ascat_total_file():
#     """Fixture to provide pyclone file for testing."""
#     return os.path.join(TEST_DATA_DIR, "ASCAT.scprofile.txt")

# @pytest.fixture
# def ascat_as_file():
#     """Fixture to provide ASCATsc allele-specific file for testing."""
#     return os.path.join(TEST_DATA_DIR, "AS_ascat_profile.txt")


# @pytest.fixture
# def read_count_file():
#     """Fixture to provide pyclone file for testing."""
#     return os.path.join(TEST_DATA_DIR, "read_counts.csv")


# @pytest.fixture
# def sample_metadata_file():
#     """Fixture to provide pyclone file for testing."""
#     return os.path.join(TEST_DATA_DIR, "sample_metadata.csv")


# @pytest.fixture
# def maf_file():
#     """Provides the path to the maf file for testing."""
#     return os.path.join(TEST_DATA_DIR, "test.maf")
