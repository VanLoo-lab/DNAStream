# DNAStream
DNAStream is a multisample, multiplatform DNA sequencing HDF5 data structure for integrated downstream evolutionary analysis. 

For the source code, visit [https://github.mdanderson.org/llweber/DNAStream](https://github.mdanderson.org/llweber/DNAStream).

See [https://pages.github.mdanderson.org/llweber/DNAStream/](https://pages.github.mdanderson.org/llweber/DNAStream/) for the API.


See [https://pages.github.mdanderson.org/llweber/DNAStream/](https://pages.github.mdanderson.org/llweber/DNAStream/) for the API.

![overview](overview.png)

## Dependencies
 - `h5py`
 - `numpy`
 - `pandas`


## Schema (under development)
The schema is currently underdevelopment is subject to change but below is the currently implemented or prospective design (*). 
```
/
 ├── SNV/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                #dataframe structure containing quality scores, number of callers, etc
     |── cluster
     |── index_map             #json string for fast loading and saving
     |-- log
 ├── sample/                     # Shared SNV index
 │   ├── labels               #short name chr:pos:ref:alt
     |── data                 #dataframe structure containing bam file path, sample code
     |── cluster
     |── index_map
     |-- log
 |-- trees/
 |   |-- SNV_trees/  
 |   |     |- trees      #edge lists of clusters
 |   |     |- data      #holds the likelihood, rank, method used to generate, etc
 |   |     |- index_map  #hold the index map
 |   |-- CNA_trees     
 |   |     |- trees      #edge lists (*) probably changed to Newick strings
 |   |     |- data      #holds the label, likelihood, rank, method used to generate, etc
 |   |     |- index_map
 |   |-- clonal_trees (joint CNA SNV tree)/
 |   |     |-- tree 
 |   |     |-- data      #holds the likelihood, rank, method used to generate, etc
 |   |     |-- index_map
 |   |     |-- genotypes (*) (structured array of node/snv/x/y/x_bar/y_bar)
 |   |     |-- clonal proportions (*) (U)
 |   |     |-- sample assignment (*)
 |   |    
 |-- copy_number/ (*)
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
 │   ├── sample_info  (*)            # Sample IDs
 │   ├── processing_parameters (*)

```