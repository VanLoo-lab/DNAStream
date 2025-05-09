import numpy as np
import h5py
from dataclasses import dataclass



#adding sample by snv metadata table

# CLINICAL_DTYPE 

@dataclass
class SampleMetadata:
    """
    Metadata schema for a sequencing sample.

    Attributes
    ----------
    sample_name : str
        Unique identifier for the sample (e.g., sample barcode or name).
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
    sample_name: str
    organism: str
    library_strategy: str
    library_source: str
    library_selection: str
    library_layout: str
    platform: str
    model: str
    center_name: str
    run: str
    study: str
    coverage: float
    modality: str
    location: str
    bam_file_path: str
    batch_id: str
    reference_build: str
    date_of_sequencing: str

    @staticmethod
    def get_dtype():
        return np.dtype([
            ("sample_name", "S100"),
            ("organism", "S40"),
            ("library_strategy", "S40"),
            ("library_source", "S40"),
            ("library_selection", "S40"),
            ("library_layout", "S10"),
            ("platform", "S40"),
            ("model", "S40"),
            ("center_name", "S60"),
            ("run", "S20"),
            ("study", "S20"),
            ("coverage", "f4"),
            ("modality", "S10"),
            ("location", "S60"),
            ("bam_file_path", h5py.string_dtype(encoding='utf-8')),
            ("batch_id", "S40"),
            ("reference_build", "S20"),
            ("date_of_sequencing", "S20"),
        ])

@dataclass
class VariantMetadata:
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
    
    chrom: str
    start_pos: int
    end_pos: int
    ref_allele: str
    alt_allele: str
    hugo: str
    entrez_gene_id: str
    variant_classification: str
    variant_type: str
    dbsnp_id: str
    filter: str
    info: str

    @staticmethod
    def get_dtype():
        return np.dtype([
            ("chrom", "S5"),
            ("start_pos", "i8"),
            ("end_pos", "i8"),
            ("ref_allele", "S10"),
            ("alt_allele", "S10"),
            ("hugo", "S40"),
            ("entrez_gene_id", "S40"),
            ("variant_classification", "S25"),
            ("variant_type", "S10"),
            ("dbsnp_id", "S20"),
            ("filter", "S20"),
            ("info", h5py.string_dtype(encoding='utf-8')),
        ])






    #("caller")   but could vary by sample
        # ("strand", "S1"),  #currently all variants report the positive strand, so removing.
            # ("context", "S10"),  #what was this?


#annotation


# SNP_DTYPE = np.dtype([
#     ("chrom", "S10"),
#     ("position", "i8"),
#     ("ref_allele", "S10"),
#     ("alt_allele", "S10"),
#     ("dbsnp_id", "S20"),
#     ("strand", "S1"),  # + or -
# ])


    #    ("pipeline_version", "S20"),