# DNAStream
DNAStream is a multisample, multiplatform DNA sequencing HDF5 data structure for integrated downstream evolutionary analysis. 

<!-- For the source code, visit [https://github.mdanderson.org/llweber/DNAStream](https://github.mdanderson.org/llweber/DNAStream). -->

See [https://pages.github.mdanderson.org/llweber/DNAStream/](https://pages.github.mdanderson.org/llweber/DNAStream/) for the API.


![overview](overview.png)


## Table of Contents
- [Dependencies](#dependencies)
- [Installation](#installation)
- [API](#build-and-host-api-locally)
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
conda create -n dnastream python=3.11
conda activate dnastream
```

Install the `dnastream` package and dependencies via pip.  

```bash
pip install .  
```

Verify the installation. 

```
python -c "from dnastream import DNAStream"
```

 Packaage is ready to use if no errors occurred!



 ## Build and host API locally
 To build and view the API locally, run the following in a terminal with the `dnastream` environment activated:
 ```bash
mkdocs build
mkdocs serve
 ```

You will see a message like `INFO    -  [10:31:55] Browser connected: http://127.0.0.1:8000/`

The click on the link or open a browser at the above address to view the API.

 ## Tutorial


### Initializing DNAStream
 Import and initialize a DNAStream object by specifying an HDF5 filename. If the file does not exist, empty datasets
 are created from the [schema](#schema-under-development). If the file exists, a connection is established to the file 
 for streaming. 
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
maf_file1 = "path/to/favorite/maf/file1.maf"
ds.add_maf_file(maf_file1)
maf_file2 = "path/to/favorite/maf/file2.maf"
#add a list of maf files
my_maf_files = [maf_file1, maf_file1]
ds.add_maf_files(my_maf_files)

 ```
 DNAStream won't add duplicate SNVs (chr:pos:ref:alt) to the index although SNV metadata will be updated at the existing indices.

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
 │   ├── labels               # Short name chr:pos:ref:alt
 │   ├── data                 # Structured array: quality scores, active
 │   ├── cluster              # Integer cluster assignments
 │   ├── log                  # Index modification log
 ├── sample/                  # Sample index
 │   ├── labels               # Sample names
 │   ├── data                 # Structured array: patient ID, source, location, file paths
 │   ├── cluster              # Integer cluster assignments
 │   ├── log                  # Index modification log
 ├── tree/
 │   ├── SNV_trees/  
 │   │   ├── trees            # Variable-length edge lists of clusters
 │   │   ├── data             # Likelihood, rank, method used to generate, etc.
 │   │   ├── labels           # Tree labels
 │   ├── CNA_trees/     
 │   │   ├── trees            # Variable-length edge lists (or Newick strings)
 │   │   ├── data             # Likelihood, rank, method used to generate, etc.
 │   │   ├── labels           # Tree labels
 ├── copy_numbers/
 │   ├── bulk/
 │   │   ├── labels           # Segment labels (chrom, start, end)
 │   │   ├── index            # Bulk-specific segment index
 │   │   ├── profile          # 3D array: (segment, sample, allele-specific CN)
 │   │   ├── logr             # 2D array: (segment, sample) logR values
 │   │   ├── baf              # 2D array: (segment, sample) B-allele frequency
 │   │   ├── log              # Modification log
 │   ├── lcm/
 │   │   ├── labels           # Segment labels
 │   │   ├── index            # LCM-specific segment index
 │   │   ├── profile          # 3D array: (segment, sample, allele-specific CN)
 │   │   ├── logr             # 2D array: (segment, sample) logR values
 │   │   ├── baf              # 2D array: (segment, sample) B-allele frequency
 │   │   ├── log              # Modification log
 │   ├── scdna/
 │   │   ├── labels           # Segment labels
 │   │   ├── index            # Single-cell segment index
 │   │   ├── profile          # (sample, segment) → allele CN tuple
 │   │   ├── log              # Modification log
 ├── read_counts/
 │   ├── variant              # 2D array: (SNV, sample) variant read counts
 │   ├── total                # 2D array: (SNV, sample) total read counts
 │   ├── log                  # Read count modifications log
 ├── metadata/
 │   ├── log                  # Metadata modifications log
 │   ├── sample_info          # Sample metadata
 │   ├── processing_parameters # Processing parameters used in analysis

```