# DNAStream
DNAStream is a multisample, multiplatform DNA sequencing HDF5 data structure for integrated downstream evolutionary analysis. 

<!-- For the source code, visit [https://github.mdanderson.org/llweber/DNAStream](https://github.mdanderson.org/llweber/DNAStream). -->

See [https://pages.github.mdanderson.org/llweber/DNAStream/](https://pages.github.mdanderson.org/llweber/DNAStream/) for the API.


![overview](overview.png)


## Table of Contents
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Tutorial](#tutorial)
  - [Initializing DNAStream](#initializing-dnastream)
  - [Adding SNVs from MAF Files](#adding-snvs-from-maf-files)
  - [Adding Read Counts](#adding-read-counts)
  - [Tree File Format for SNV Phylogenies](#tree-file-format-for-snv-phylogenies)
  - [Adding SNV Phylogenies](#adding-snv-phylogenies)
  - [Viewing Modification Logs](#viewing-modification-logs)
- [Schema (Under Development)](#schema-under-development)


## Dependencies
 - `h5py`
 - `numpy`
 - `pandas`


## Installation

Clone the git repository and create a conda environment (recommended).  

```bash
git clone git@github.mdanderson.org:llweber/DNAStream.git
cd DNAstream
mamba create -n dnastream python=3.11
mamba activate dnastream
```

Install the `dnastream` package and dependencies via pip.  

```bash
pip install .  
```

Verify the installation. 

```
python -c "from dnastream import DNAStream"
```

 Packing is ready to use if no errors occurred. 


 ## Tutorial


### Initializing DNAStream
 Import and initialize a DNAStream.
 ```python
#create a DNAStream object and corresponding HFD5 file
from dnastream import DNAStream

#use verbose mode to see additional processing messages
ds = DNAStream("myfile.h5", verbose=True)

```

### Adding SNVs from MAF Files
Add SNVs to index and associated metadata from a MAF file(s).
```python
#add SNVs to index and associated metadata with 
maf_file = "path/to/favorite/maf/file"
ds.add_maf_file(maf_file)

#add a list of maf files
my_maf_files = [maf_file, maf_file]
ds.add_maf_file(my_maf_files)
#DNAStream won't add duplicate SNVs (chr:pos:ref:alt) to the index. 
 ```

### Adding Read Counts
Add read counts for single-cell data
```python
#must have columns ordered snv | sample | variant | total
#column names will be ignored but the order matters. 
read_count_file = "my_read_count_file.csv"


ds.add_read_counts(read_count_file, source="scdna")

#any new SNVs or samples not in the index will be added 
#retreive the index logs as Pandas DataFrames to check index changes
snv_log = ds.get_snv_log()
print(snv_log)
sample_log = ds.get_sample_log()
print(sample_log)
```

### Tree File Format for SNV Phylogenies
Example of tree file format for SNV phylogenies:
```
#tree 1
8	4
4	6
4	0
6	3
6	7
8	2
8	5
#tree 2
8	4
4	6
4	0
6	7
8	2
8	5
4	3
#tree 3
8	4
4	6
4	0
6	3
8	2
8	5
4	7
#tree 4
8	4
4	6
4	0
8	2
8	5
4	3
4	7

```

### Add SNV phylogenies
Add SNV phylogenies to DNAStream.
```python
#ensure file is in the above format
tree_file = "path/to/trees.txt"
method= "conipher"
ds.add_trees_from_file(tree_file, tree_type="SNV", method=method, safe=True)
```

Note: Safe mode checks to see if trees have already been added from the source
file.  If safe=True, trees will be not be added.  If safe=False, duplicate
trees may be added if trees were previously loaded from the same source file. 


### Viewing Modification Logs
View the modification log that tracks changes to any dataset.
```python
log = ds.get_dataset_log()
print(log)
```


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