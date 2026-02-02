# DNAStream Core Object Model

## Entities
Entities are first-class objects with stable, immutable identifiers stored in registries.  
1.	**Sample**: A biological or computational unit of observation (e.g., bulk sample, single cell, pseudobulk, LCM spot, ctDNA).  
2.	**Variant**: A genomic variant, e.g., SNV, indel, defined by genomic coordinates and alleles.  
3.	**Segmentation**: A specific genomic segmentation produced for a given modality and analysis. eg., bulk CN segmentation, scDNA segmentation.  
4.	**Segment**: A single genomic interval belonging to a specific Segmentation.  
5.	**Bin**: An optional finer-grained genomic interval belonging to a Segmentation, which may be aggregated into Segments.   
6.	**Result**: An immutable output of an analysis, e.g., SNV clusterings, clone trees, copy number calls, mixture proportion estimates, that references one or more entities and/or measurements with explicit provenance and versioning.  

## Data categories
These are stored data objects keyed by entities, not registries themselves.
1. **Measurements**: Observed or directly computed data derived from sequencing (e.g., read counts, BAFs, logR), stored as typed tensors indexed by entities. Measurements cannot be designated canonical independently of the entities they index.

## Operations
1.	**register**: add/read registries with stable IDs
2.	**attach**: store measurement tensors keyed by registries
3.	**link**: store mappings between registries, e.g. SNV to segment, bin to segment
4.	**publish**: create immutable results with run metadata
5.	**designate**: set canonical references as pointers