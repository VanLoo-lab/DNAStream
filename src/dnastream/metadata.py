import numpy as np
import h5py
from .datatypes import STR_DTYPE

SAMPLE_DTYPE = np.dtype(
    [
        ("sample_name", "S100"),  #unique identifier for the sample
        ("organism", "S40"),  #e.g. human, mouse
        ("library_strategy", "S40"),  #e.g.  genomic
        ("library_source", "S40"),  #type of source material being sequenced e.g. DNA
        ("library_selection", "S40"),  #method used to select/enrich the library
        ("library_layout", "S10"),  #paired or single end
        ("platform", "S40"), #type of sequencing platform
        ("model", "S40"),   #sequencing model
        ("center_name", "S60"),  #name of sequencing center
        ("run", "S20"),
        ("study", "S20"),
        ("coverage", "f4"),   #mean sequence coverage
        ("modality", "S10"),  #bulk, scdna, lcm
        ("location", "S60"),  #short id of where sample originated e.g. LR_1
        ("bam_file_path", h5py.special_dtype(vlen=str)),
        ("batch_id", "S40"),
        ("reference_build", "S20"),  #e.g. hg19, hg38
        ("date_of_sequencing", "S20"),
    ]
)


SNV_DTYPE = np.dtype([
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
    ("strand", "S1"),
    ("filter", "S20"),
    ("context", "S10"),
    ("info", h5py.special_dtype(vlen=str)),
])


SNP_DTYPE = np.dtype([
    ("chromosome", "S10"),
    ("position", "i8"),
    ("ref_allele", "S10"),
    ("alt_allele", "S10"),
    ("dbsnp_id", "S20"),
    ("strand", "S1"),
])


    #    ("pipeline_version", "S20"),