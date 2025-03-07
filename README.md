# DNAStream: A multisample, multiplatform DNA sequencing HDF5 data structure for integrated downstream evolutionary analsysis 

See [https://pages.github.mdanderson.org/llweber/DNAStream/](https://pages.github.mdanderson.org/llweber/DNAStream/) for the API.

## Dependencies
 - `h5py`
 - `numpy`
 - `pandas`



## Schema (under development)
```
/
 ├── SNV/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing quality scores, number of callers, etc
     |── cluster
     |── index_map             #json string for fast loading and saving
     |-- log
 ├── sample/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing bam file path, sample code
     |── cluster
     |── index_map
     |-- log
 ├── edge_list/               # Read count matrices
 │   ├── bulk/                  # Bulk sequencing read counts
 │   │   ├── variant       # SNVs x Samples (variant read counts)
 │   │   ├── total         # SNVs x Samples (total read counts)
 │   ├── lcm/                    # LCM sequencing read counts
 │   │   ├── variant       
 │   │   ├── total         
 │   ├── scdna/                  # scDNA-seq read counts
 │   │   ├── variant       
 │   │   ├── total 
 |-- trees/
 |   |-- SNV_trees/   #edge lists of clusters
 |   |     |- trees
 |   |     |- data      #holds the likelihood, rank, method used to generate, etc
 |   |     |- index_map  #hold the index map
 |   |-- CNA_trees     # cell lineage trees   #newick strings
 |   |     |- trees
 |   |     |- data      #holds the likelihood, rank, method used to generate, etc
 |   |     |- index_map
 |   |-- clonal_trees (joint CNA SNV tree)/
 |   |     |-- tree (edge list)
 |   |     |-- genotypes  (structured array of node/snv/x/y/x_bar/y_bar)
 |   |     |-- clonal proportions (U)
 |   |     |-- sample assignment 
 |   |    
 |-- copy_number/
 |   ├── /bulk/
 |   │   ├── /segments    # Bulk-specific segment index
 |   │   ├── /profiles       # Tensor: (sample, segment, allele CN, proportion μ)
 |   │   ├── /metadata    # Bulk-specific metadata
 |   |   |-- /log          # logging
 |   │
 |   ├── /single_cell/
 |   │   ├── /segments    # Single-cell segment index
 |   │   ├── /profiles       # (sample, segment) → allele CN tuple
 |   │   ├── /metadata    # Single-cell metadata
 |   │
 |   ├── /lcm/
 |   │   ├── /segments    # LCM-specific segment index
 |   │   ├── /profiles       # (sample, segment) → allele CN tuple
 |   |   |--/logR
 |   |   |--/baf
 |   │   ├── /metadata    # LCM-specific metadata
 ├── metadata/ 
     |-- log                  # Metadata storage
 │   ├── sample_info              # Sample IDs
 │   ├── processing_parameters

```