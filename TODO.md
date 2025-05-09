# TODO list

## Debugging
- FIXED (05/08/2025): Duplicate SNVs are being added to the index when files are loaded multiple times.  
   The issue is that the _index_cache is not loaded, so duplicates are always added because the cache is not intialized from the stored metadata.
   It turns out labels name were stored as bytes and the new labels were stored in unicode and therefore were registering as being in the label_cache or file.
   The fix was to decode the stored labels from bytes to unicode when loading the index cache.

## Testing  
<!-- - add test script for index_manager.py -->
- test adding pseudobulk layers
- generate test maf file
<!-- - test adding a global index -->


## I/O Functions
<!-- - parse_battenberg_file -->
<!-- - parse_pyclone -->
- parse_vcf_file
- add column mappings for MAF 
<!-- - parse ascat output -->

## Data structure
- add layer for lcm coordinates
<!-- - store path to bam files for each sample -->
- fix metadata management
<!-- - modify constants to use Enums in DNAStream -->



<!-- - map global indices to controlled datasets and have index changes trigger resizing of controlled tables -->


## API/Readme
- Add documentation for data structure to API
- add example for adding SNV trees by edge lists 
- add example for adding a global index (only allow new tables to be controlled for now.)
- fix docs for factory functions accessing GlobalIndex

## Indices 
   - Behavior if an SNV appears in multiple MAF files.  Should there be a data table with SNV by sample that holds sample specific metadata 
   such as base quality scores...


## Develop scripts for loading files into DNAStream object