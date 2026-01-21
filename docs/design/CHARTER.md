# DNAStream Charter

DNAStream is an effiecient storage + access layer for processed tumor DNA sequencing data used in intra-tumor heterogeneity and tumor phylogenetic analysis, with reproducible provenance and safe collaboration. The goal of DNAStream is develop a reliable, scalable, and maintainable data systems for cancer genomics data.

## What DNAStream is not
- DNAStream is not a storage container for raw data, such as fastq and BAM files
- DNAStream is not a workflow pipeline like Snakemake or Nextflow
- DNAStream is not a visualization tool
- DNAStream is not a full analysis tool but does store the input and output of relevant anlyses

## DNAStream design requirements
1.	**Multi-user safe access** that supports concurrent readers and controlled write access without data corruption (e.g., single-writer with append-only artifacts).
2.	**Explicit immutability rules**: *core entities*, e.g., samples, SNVs, segmentations, and *analysis artifacts*, e.g., trees,copy numbers, model fits, are immutable; updates create new versions or views rather than overwriting existing data.
3.	**Explicit versioning of schema and analysis outputs** to support reproducibility, comparison across runs, and backward compatibility.
4.	**Identity-based indexing and access**: all key entities use stable, immutable identifiers, and data access must never rely solely on positional meaning.
5.	**Identity-preserving subsetting**: subsetting or reordering data preserves entity identifiers and referential integrity across all associated data structures.
6.	**Separation of concerns**: measured data (e.g., read counts, BAFs), derived analyses (e.g., copy number states, cluster assignments, trees), and metadata are stored distinctly and referenced explicitly.
7.	**Canonical reference designation**: DNAStream supports explicit designation of canonical reference entities (e.g., SNV sets, sample sets, CN profiles, trees) with versioning, provenance, and deprecation rather than silent replacement.
8.	**Efficient data storage and access** with chunked, lazy reading suitable for large-scale genomic data.
9.	**Efficient I/O interoperability** tailored to the outputs of common bioinformatic tools and pipelines (e.g., PyClone, Battenberg, ASCAT).
10.	**Tidy metadata conventions** with explicit, well-typed fields suitable for programmatic access and filtering.
11.	**Schema versioning** that precisely specifies the underlying data structure and governs compatibility and migration.
12.	**Validation on write and modification** to prevent violations of schema, identity, or design invariants.
13.	**Fail loudly and early**: invalid shapes, mismatched indices, missing identifiers, or schema violations raise explicit errors rather than being silently coerced.



## Use cases
1. Load SNV read counts for a subset of samples, compute VAF; filter SNVs, write filtered view without rewriting raw.
2. Store multiple CN segmentations (bulk/LCM/scDNA); query major/minor CN for a sample and region.
3. Link SNVs and segments for a given segmentation and export joined tables for downstream modeling.
4. Store candidate trees with scores; reproduce ranking later with provenance.
5. Attach analysis outputs (clusters, mixture proportions) with clear ownership and run metadata.
6. Two analysts can add new results without clobbering each other (single-writer lock + append-only artifacts).
7. Ability to easily load output data from common bioinformatics tools and pipelines 
8. Set the source of ground truth for SNVs, quality samples in analysis and copy number profiles so all analysts are using the most up to date data.
