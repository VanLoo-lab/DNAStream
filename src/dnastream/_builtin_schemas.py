import re
import h5py
import numpy as np
from .schema import Schema, Field
from .constants import EVENTS, STR_DTYPE, SCHEMA_VERSION


def str_validator(x):
    pass


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


SAMPLE_REGISTRY_FIELDS = REGISTRY_SPINE + (
    Field(
        "sample_name",
        STR_DTYPE,
        True,
        str_validator,
    ),
    Field("organism", STR_DTYPE, str_validator),
    Field("library_strategy", STR_DTYPE, str_validator),
    Field("library_source", STR_DTYPE, str_validator),
    Field("library_selection", STR_DTYPE, str_validator),
    Field("library_layout", "S10"),
    Field("platform", STR_DTYPE, str_validator),
    Field("model", STR_DTYPE, str_validator),
    Field("center_name", STR_DTYPE, str_validator),
    Field("run", STR_DTYPE, str_validator),
    Field("study", STR_DTYPE, str_validator),
    Field("coverage", "f4"),
    Field("modality", STR_DTYPE, str_validator),
    Field("location", STR_DTYPE, str_validator),
    Field("bam_file_path", STR_DTYPE, str_validator),
    Field("batch_id", STR_DTYPE, str_validator),
    Field("reference_build", STR_DTYPE, str_validator),
    Field("date_of_sequencing", STR_DTYPE, str_validator),
)


def builder_sample_label(sample_name: str) -> str:
    s = str(sample_name).strip()
    s = re.sub(r"\s+", " ", s)
    return s


SAMPLE_SCHEMA = Schema(
    fields=SAMPLE_REGISTRY_FIELDS,
    version=SCHEMA_VERSION,
    label_from=("sample_name",),
    label_required=True,
    label_builder=builder_sample_label,
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
VARIANT_REGISTRY_FIELDS = REGISTRY_SPINE + (
    Field("chrom", STR_DTYPE, False, str_validator),
    Field("start_pos", "i8", False),
    Field("end_pos", "i8", False),
    Field("ref_allele", STR_DTYPE, False, str_validator),
    Field("alt_allele", STR_DTYPE, False, str_validator),
    Field("hugo", STR_DTYPE, False, str_validator),
    Field("entrez_gene_id", STR_DTYPE, False, str_validator),
    Field("variant_classification", STR_DTYPE, False, str_validator),
    Field("variant_type", STR_DTYPE, False, str_validator),
    Field("dbsnp_id", STR_DTYPE, False, str_validator),
    Field("filter", STR_DTYPE, False, str_validator),
    Field("info", STR_DTYPE, False, str_validator),
    Field("source", STR_DTYPE, False, str_validator),  # method
)


def normalize_chrom(chrom: str) -> str:
    c = str(chrom).strip()
    c = c[3:] if c.lower().startswith("chr") else c
    c = c.upper()
    if c == "MT":
        c = "M"
    return f"chr{c}"


def build_variant_label(chrom, start_pos, ref_allele, alt_allele) -> str:
    c = normalize_chrom(chrom)
    pos = int(start_pos)
    ref = str(ref_allele).strip().upper()
    alt = str(alt_allele).strip().upper()
    if not ref or not alt:
        raise ValueError("ref_allele and alt_allele must be non-empty")
    return f"{c}:{pos}:{ref}:{alt}"


VARIANT_SCHEMA = Schema(
    VARIANT_REGISTRY_FIELDS,
    version=SCHEMA_VERSION,
    label_from=("chrom", "start_pos", "ref_allele", "alt_allele"),
    label_builder=build_variant_label,
    label_required=True,
    label_normalizer=None,
)


SNP_FIELDS = REGISTRY_SPINE + (
    Field("chrom", "S10", True, STR_DTYPE),
    Field("start_pos", "i8", True, STR_DTYPE),
    Field("ref_allele", "S10", True, STR_DTYPE),
    Field("alt_allele", "S10", True, STR_DTYPE),
    Field("dbsnp_id", "S20", False, STR_DTYPE),
    Field("strand", "S1", False, None),  # + or -
)

SNP_SCHEMA = Schema(
    SNP_FIELDS,
    version=SCHEMA_VERSION,
    label_from=("chrom", "start_pos", "ref_allele", "alt_allele"),
    label_builder=build_variant_label,
    label_required=True,
    label_normalizer=None,
)

# lambda x: raise ValueError("strand must be either + or -") if x not in {"+", "-"}


REGISTRY_SCHEMAS = {
    "sample": SAMPLE_SCHEMA,
    "variant": VARIANT_SCHEMA,
    "snp": SNP_SCHEMA,
}


# annotation


#    ("pipeline_version", "S20"),


def validate_event(x):
    if x not in EVENTS:
        raise ValueError(f"Event '{x}, not valid must be ones of {','.join(EVENTS)}")


PROVENANCE_LOG = (
    Field("id", STR_DTYPE, True, None),
    Field("timestamp", STR_DTYPE, True, None),  # ISO8601 Z
    Field("user", STR_DTYPE, True, None),
    Field("scope", STR_DTYPE, True, None),
    Field("event", STR_DTYPE, True, validate_event),
    Field("dataset", STR_DTYPE, True, None),  # full HDF5 path
    Field("source", STR_DTYPE, False, None),  # module.qualname (optional)
    Field("info", STR_DTYPE, False, None),  # JSON string (optional)
)

PROVENANCE_LOG_SCHEMA = Schema(
    version=SCHEMA_VERSION, fields=PROVENANCE_LOG, label_required=False
)

PROVENANCE_SCHEMAS = {"log": PROVENANCE_LOG_SCHEMA}
