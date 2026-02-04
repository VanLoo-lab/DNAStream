# DNAStream

DNAStream is an HDF5-backed, multi-modal data structure for organizing DNA sequencing data and downstream evolutionary analysis. It provides compact on-disk storage, fast partial reads, and a structured way to track entities, links, and changes over time.

## Beta scope

This beta release focuses on:
- **Registry**: typed registries with built-in schemas for core entities (e.g., samples, variants, SNPs) and activation status
- **Provenance**: lightweight event logging for dataset changes (**create**, **append**, **modify**, **designate**)

Expect the API and on-disk layout to evolve during beta.

![overview](assets/dnastream.png)
*Beta includes Registry + Provenance. Measurements and Results are planned (marked with * in the diagram).*

## Key features

- **Efficient storage and access** via a chunked HDF5 file with lazy reads for large cohorts
- **Entity registries with schemas** to validate fields, manage activation, and support consistent linking across datasets and analyses
- **Provenance logging** of change events to support reproducibility and collaboration

## Coming soon

- **Measurements** linked to registered entities (e.g., variant/total read counts, binned counts)
- **Results** storage and retrieval (e.g., copy number calling, clonal trees)
- **Canonical result pointers** to mark the active/preferred outputs among multiple runs
- **Custom schemas** for specialized registries, measurements, and results
- **Multi-user workflows** with a clear concurrency policy for write access



<!-- ## Table of Contents
- [Dependencies](#dependencies)
  - [Optional dependences](#optional-dependencies)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [Documentation](#documentation)
- [Unit tests](#unit-tests) -->






