from enum import Enum


class GlobalIndexName(Enum):
    SNV = "SNV"
    SAMPLE = "sample"
    SNP = "SNP"


class LocalIndexName(Enum):
    TREES_SNV = "trees_SNV"
    TREES_CNA = "trees_CNA"
    TREES_CLONAL = "trees_CLONAL"
    COPY_NUMBERS_BULK = "copy_numbers_bulk"
    COPY_NUMBERS_LCM = "copy_numbers_lcm"
    COPY_NUMBERS_SCDNA = "copy_numbers_scdna"


class Modalities(Enum):
    BULK = "bulk"
    LCM = "lcm"
    SCDNA = "scdna"


class TreeType(Enum):
    SNV = "SNV"
    CNA = "CNA"
    CLONAL = "CLONAL"


class SchemaGroups(Enum):
    TREES = "trees"
    COPY_NUMBERS = "copy_numbers"
    READ_COUNTS = "read_counts"
    ALLELE_COUNTS = "allele_counts"
