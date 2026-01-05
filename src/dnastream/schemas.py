import re
import h5py
import numpy as np
import hashlib
import json

SCHEMA_VERSION = "0.1.0"
"""
Metadata schema for a sequencing sample.

Attributes
----------
organism : str
    Organism from which the sample was derived (e.g., 'human', 'mouse').
library_strategy : str
    Sequencing strategy used (e.g., 'WGS', 'WES', 'RNA-Seq', 'genomic').
library_source : str
    Type of source material being sequenced (e.g., 'genomic', 'transcriptomic', 'metagenomic', 'DNA', 'RNA').
library_selection : str
    Method used to enrich or select the library (e.g., 'PCR', 'hybrid capture', 'polyA', 'random').
library_layout : str
    Sequencing layout (e.g., 'paired-end', 'single-end').
platform : str
    Type of sequencing platform (e.g., 'Illumina', 'PacBio', 'Nanopore').
model : str
    Specific model of the sequencing instrument (e.g., 'NovaSeq6000', 'NextSeq', 'HiSeq4000').
center_name : str
    Name of the sequencing center or facility where sequencing was performed.
run : str
    Sequencing run identifier (e.g., run accession or batch number).
study : str
    Study or project identifier (e.g., project code or accession).
coverage : float
    Mean sequence coverage (depth) of the sample (e.g., 30.0 for 30x coverage).
modality : str
    Sequencing modality (e.g., 'bulk', 'scdna', 'lcm', 'single-cell').
location : str
    Sample origin or location (e.g., tissue, slide position, well, or coordinate).
bam_file_path : str
    Path to the BAM file corresponding to the sample.
batch_id : str
    Identifier for the experimental or sequencing batch (e.g., batch number or date).
reference_build : str
    Reference genome build used for alignment (e.g., 'hg19', 'hg38', 'mm10').
date_of_sequencing : str
    Date the sample was sequenced (e.g., '2023-06-15').
"""

STR_DTYPE = h5py.string_dtype("utf-8")

REGISTRY_SPINE = (
    ("id", STR_DTYPE),  # UUIDv4 string
    ("label", STR_DTYPE),  # user-facing unique key
    ("idx", np.int64),
    ("active", np.bool_),
    ("created_at", STR_DTYPE),  # ISO8601 Z
    ("created_by", STR_DTYPE),
)

SAMPLE_REGISTRY_DTYPE_SPEC = REGISTRY_SPINE + (
    ("sample_name", STR_DTYPE),
    ("organism", STR_DTYPE),
    ("library_strategy", STR_DTYPE),
    ("library_source", STR_DTYPE),
    ("library_selection", STR_DTYPE),
    ("library_layout", "S10"),
    ("platform", STR_DTYPE),
    ("model", STR_DTYPE),
    ("center_name", STR_DTYPE),
    ("run", STR_DTYPE),
    ("study", STR_DTYPE),
    ("coverage", "f4"),
    ("modality", STR_DTYPE),
    ("location", STR_DTYPE),
    ("bam_file_path", STR_DTYPE),
    ("batch_id", STR_DTYPE),
    ("reference_build", STR_DTYPE),
    ("date_of_sequencing", STR_DTYPE),
)


"""
Metadata schema for a single nucleotide variants (SNV) and single nucleotide polymorphisms (SNP).

Attributes
----------
chrom : str
    Chromosome identifier (e.g., 'chr1').
pos : int
    Genomic start position of the SNV.
end_pos : int
    Genomic end position of the SNV.
ref_allele : str
    Reference allele.
alt_allele : str
    Alternate allele.
hugo : str
    HUGO gene symbol.
gene : str
    Entrez gene ID.
filter : str
    Variant filter status (e.g., 'PASS').
variant_classification : str
    Functional effect (e.g., 'Missense_Mutation').
variant_type : str
    Variant type (e.g., 'SNP', 'INS').
dbsnp_id : str
    dbSNP reference ID.
info : str
    Additional metadata (e.g., from INFO field in VCF/MAF).
"""
VARIANT_REGISTRY_DTYPE_SPEC = REGISTRY_SPINE + (
    ("chrom", STR_DTYPE),
    ("start_pos", "i8"),
    ("end_pos", "i8"),
    ("ref_allele", STR_DTYPE),
    ("alt_allele", STR_DTYPE),
    ("hugo", STR_DTYPE),
    ("entrez_gene_id", STR_DTYPE),
    ("variant_classification", STR_DTYPE),
    ("variant_type", "S10"),
    ("dbsnp_id", STR_DTYPE),
    ("filter", STR_DTYPE),
    ("info", STR_DTYPE),
)


def dtype_from(spec):
    return np.dtype(list(spec))


def columns_from(spec):
    return [k for k, _ in spec]


def schema_hash(spec):
    payload = json.dumps(
        [(k, str(v)) for k, v in spec],
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def compile_schema(
    spec,
    *,
    label_from=None,
    label_normalizer=None,
    label_required: bool = False,
):
    spec = tuple(spec)

    if label_from is not None:
        if isinstance(label_from, (list, tuple)):
            label_from = tuple(label_from)
        else:
            raise TypeError("label_from must be a tuple/list of field names or None")

    return {
        "spec": spec,
        "dtype": dtype_from(spec),
        "columns": columns_from(spec),
        "schema_pairs": [(k, str(v)) for k, v in spec],
        "schema_hash": schema_hash(spec),
        "schema_version": SCHEMA_VERSION,
        "label_from": label_from,
        "label_normalizer": label_normalizer,
        "label_required": bool(label_required),
    }


def normalize_sample_label(sample_name: str) -> str:
    s = str(sample_name).strip()
    s = re.sub(r"\s+", " ", s)
    return s


def normalize_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    c = c[3:] if c.lower().startswith("chr") else c
    c = c.upper()
    if c == "MT":
        c = "M"
    return f"chr{c}"


def normalize_variant_label(chrom, start_pos, ref_allele, alt_allele) -> str:
    c = normalize_chrom(chrom)
    pos = int(start_pos)
    ref = str(ref_allele).strip().upper()
    alt = str(alt_allele).strip().upper()
    if not ref or not alt:
        raise ValueError("ref_allele and alt_allele must be non-empty")
    return f"{c}:{pos}:{ref}:{alt}"


SAMPLE_REGISTRY = compile_schema(
    SAMPLE_REGISTRY_DTYPE_SPEC,
    label_from=("sample_name",),  # fields used to compute label
    label_normalizer=normalize_sample_label,
    label_required=True,
)

VARIANT_REGISTRY = compile_schema(
    VARIANT_REGISTRY_DTYPE_SPEC,
    label_from=("chrom", "start_pos", "ref_allele", "alt_allele"),
    label_normalizer=normalize_variant_label,
    label_required=True,
)


REGISTRIES = {"sample": SAMPLE_REGISTRY, "variant": VARIANT_REGISTRY}


# annotation


# SNP_DTYPE = np.dtype([
#     ("chrom", "S10"),
#     ("position", "i8"),
#     ("ref_allele", "S10"),
#     ("alt_allele", "S10"),
#     ("dbsnp_id", "S20"),
#     ("strand", "S1"),  # + or -
# ])


#    ("pipeline_version", "S20"),
