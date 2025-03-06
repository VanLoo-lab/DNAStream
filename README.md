# DNAStream: A multisample, multiplatform DNA sequencing HDF5 data structure for integrated downstream evolutionary analsysis 

## Dependencies
 - `h5py`
 - `numpy`
 - `pandas`

## Hierarchy
 ├── SNV/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |-- data                 #dataframe structure containing quality scores, number of callers, etc
     |-- cluster
 ├── sample/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |-- data                 #dataframe structure containing bam file path, sample code
     |-- cluster
 ├── read_counts/               # Read count matrices
 │   ├── bulk/                  # Bulk sequencing read counts
 │   │   ├── variant_reads       # SNVs x Samples (variant read counts)
 │   │   ├── total_reads         # SNVs x Samples (total read counts)
 │   ├── lcm/                    # LCM sequencing read counts
 │   │   ├── variant_reads       
 │   │   ├── total_reads         
 │   ├── scDNA/                  # scDNA-seq read counts
 │   │   ├── variant_reads       
 │   │   ├── total_reads         
 ├── metadata/                   # Metadata storage
 │   ├── sample_info              # Sample IDs
 │   ├── processing_parameters