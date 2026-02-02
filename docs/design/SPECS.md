# DNAStream Specifications (Draft)

This document is the on-disk contract for DNAStream: what MUST exist in a DNAStream file, what MAY exist, and how readers/writers interpret it.

This spec is intentionally minimal at first. It defines the “spine” of the format so the project can grow without breaking invariants.

---

## 0. Scope

DNAStream stores **processed** tumor DNA sequencing data and analysis outputs for intra-tumor heterogeneity and phylogenetic analysis.

Out of scope:
- raw FASTQ/BAM/CRAM content (may be referenced by path/URI)
- workflow orchestration
- visualization

---

## 1. File identity and required header

A DNAStream file MUST be an HDF5 file with the following required file-level attributes:

- `dnastream_format` (string)
  - MUST equal: `"DNAStream"`
- `schema_version` (string)
  - Semantic version: `MAJOR.MINOR.PATCH` (e.g., `"0.1.0"`, `"1.0.0"`)
- `created_at` (string)
  - ISO-8601 UTC timestamp (e.g., `"2025-12-17T16:02:10Z"`)
- `created_by` (string)
  - free-form user identifier (e.g., OS user, email)

Recommended file-level attributes (SHOULD):
- `file_uuid` (string)
- `dnastream_software_version` (string)

### 1.1 Schema compatibility rules

- Readers MUST refuse to open files with an unsupported **major** `schema_version`.
- Minor/Patch versions MUST be backward compatible unless explicitly documented.
- Breaking changes require explicit migration. Silent automatic migration is prohibited.

---

## 2. Top-level namespaces

DNAStream reserves the following top-level groups. A compliant file MUST contain them (they may be empty unless otherwise specified).

- `/registry` — entity registries with stable IDs
- `/measurements` — measured tensors keyed by entity IDs
- `/links` — explicit mappings between registries
- `/results` — immutable published analysis outputs
- `/canonical` — pointers designating canonical entities/results
- `/provenance` — append-only run and modification records

Notes:
- This layout is the “clean” target model. If legacy layouts exist (e.g., `/SNV`, `/sample`, `/copy_numbers`), they MUST be treated as legacy schemas with explicit `schema_version` and migration.

---

## 3. Entity registries

### 3.0 Registry storage format (RESOLVED for v0.1)

For schema version v0.1, registries MUST be stored as **compound (row-wise) HDF5 datasets**.

- Each registry is a single dataset at `/registry/<entity_name>`.
- Each row corresponds to one entity.
- Each field corresponds to a column.
- Writers MUST NOT store registries as multiple per-column datasets in v0.1.

Rationale:
- Compound registries make alignment invariants (all columns same length, atomic row append) mechanically enforceable.
- Registry tables are expected to be small relative to measurement tensors; reading full rows is acceptable.

Future extension (non-normative): a later schema MAY add optional columnar mirrors for performance/interoperability, but the compound dataset remains the source of truth.

Registries are tables that define **first-class entities** and their stable, immutable identifiers.

### 3.1 Common registry requirements

Each registry MUST:
- be stored under `/registry/<entity_name>`
- provide a stable ID column named `<entity>_id` (string)
- `<entity>_id` MUST be a UUID (UUIDv4 canonical text form, lowercase hex with hyphens, e.g., `"550e8400-e29b-41d4-a716-446655440000"`).
- `label` is a human-readable identifier and MAY be non-unique; code MUST NOT rely on `label` for identity.
- Registries SHOULD provide an optional `external_id` (string) column to store tool-specific IDs (e.g., VCF IDs, PyClone mutation IDs) when needed.
- preserve IDs across subsetting and reordering
- include audit fields for provenance

Each registry MUST include these columns (minimum):
- `<entity>_id` (string)
- `label` (string)
- `active` (bool)
- `created_at` (string; ISO-8601 UTC)
- `created_by` (string)
- `run_id` (string; may be empty if not tied to a run)

Registry rows MAY include domain-specific columns.

### 3.2 Required registries (minimum viable)

A DNAStream file MUST include the following registries:

1) `/registry/samples`
- required additional columns:
  - `modality` (string; e.g., `bulk`, `lcm`, `scdna`, `ctdna`, `pseudobulk`)
  - `patient_id` (string; may be empty)

2) `/registry/variants`
- required additional columns:
  - `chrom` (string)
  - `pos` (int)
  - `ref` (string)
  - `alt` (string)

3) `/registry/segmentations` (may be empty initially)
- required additional columns:
  - `modality` (string)
  - `build` (string; e.g., `hg38`)

4) `/registry/segments` (may be empty initially)
- required additional columns:
  - `segmentation_id` (string)
  - `chrom` (string)
  - `start` (int)
  - `end` (int)

5) `/registry/results` (may be empty initially)
- required additional columns:
  - `result_type` (string)
  - `status` (string; e.g., `draft`, `published`, `deprecated`)

Optional registries:
- `/registry/bins`

### 3.3 Registry immutability (v0.1)

Registries define entity identity and therefore have strict immutability rules.

- Registry rows are append-only: writers MAY append new rows but MUST NOT modify existing rows in-place.
- Deactivation/deprecation MUST be represented as a new registry **version** (new dataset) rather than editing `active` in-place.

Registry versioning (path convention):
- For v0.1, a file MUST use exactly one of the following conventions consistently:
  1) `/registry/<entity_name>__v<semver>` (example: `/registry/variants__v0.1.0`)
  2) `/registry/<entity_name>/v<semver>` (example: `/registry/variants/v0.1.0`)

Canonical pointers SHOULD be used to designate the current registry version (e.g., `/canonical/variants_registry`).

---

## 4. Measurements

Measurements are typed arrays/tensors keyed by registry IDs. They are not registries.

### 4.1 Measurement declaration

For schema version v0.1 (Option H1), the canonical on-disk representation of a measurement is **append-by-block**.

Each measurement MUST be stored under:

- `/measurements/<name>/blocks/index` — a block index table
- `/measurements/<name>/blocks/<block_id>/data` — one or more write-once block datasets

A measurement MAY additionally provide a convenience dataset:
- `/measurements/<name>/data` — an optional fully materialized tensor (cache/view)

If both `/blocks/*` and `/data` are present, `/blocks/*` is the source of truth. Readers MAY ignore `/data`.

Each measurement group MUST include attributes describing axes:

- `axes` (list of strings)
  - Example: `["variant_id", "sample_id"]`
- `dtype_class` (string)
  - Examples: `counts`, `float`, `categorical`, `interval`

Recommended attributes:
- `units` (string)
- `description` (string)
- `source` (string; tool/pipeline name)

### 4.1.1 Measurement immutability and growth (v0.1)

Measurements are immutable once written, with one permitted growth pattern: **append-by-block** along a declared growable axis.

- In-place modification of existing measurement values is prohibited.
- If a measurement needs to change (e.g., filtering, normalization, correction), writers MUST publish a new measurement name or a new measurement version.
- Measurements MAY declare one growable axis via attribute `grow_axis` (string; one of the axis names listed in `axes`).
  - If `grow_axis` is absent, the measurement is fully write-once (no growth allowed).

Canonical layout (v0.1): append-by-block
- Base path: `/measurements/<name>/`
- Required group:
  - `/measurements/<name>/blocks/`
- Required dataset:
  - `/measurements/<name>/blocks/index` — a compound table describing each block with columns:
    - `block_id` (string; UUIDv4 recommended)
    - `start` (int)
    - `length` (int)
    - `run_id` (string)
    - `created_at` (string; ISO-8601 UTC)
    - `created_by` (string)

Block data storage (v0.1):
- Writers SHOULD store each block as a separate write-once dataset at:
  - `/measurements/<name>/blocks/<block_id>/data`
- Readers MUST reconstruct the logical tensor by concatenating blocks in ascending `start` order along `grow_axis`.

Optional `/data` dataset:
- Writers MAY store `/measurements/<name>/data` as a materialized cache for fast reads.
- If `grow_axis` is present, writers MUST NOT update `/data` in-place; instead they MUST either omit `/data` or rewrite it as a new measurement name/version.

Notes:
- Blocks themselves are write-once.
- Writers MUST ensure blocks do not overlap in `[start, start+length)` along `grow_axis`.

### 4.2 Required initial measurements (suggested)

These are not required for an empty file, but are the initial targets for interoperability:

- `snv/alt_counts` with axes `["variant_id", "sample_id"]`
- `snv/total_counts` with axes `["variant_id", "sample_id"]`

Copy number-related measurements SHOULD be keyed to `segment_id` (and therefore imply a segmentation choice):
- `cn/logr` with axes `["segment_id", "sample_id"]`
- `cn/baf` with axes `["segment_id", "sample_id"]`
- `cn/major_minor` (or allele-specific) with axes `["segment_id", "sample_id", "allele"]`

---

## 5. Links

Links are explicit mapping tables between registries.

Stored under:
- `/links/<name>`

Each link MUST declare:
- `from` (string; entity id type)
- `to` (string; entity id type)
- `schema` (string; short description)

Common links:
- variant → segment (for a given segmentation)
- bin → segment

---

## 6. Results

Results are immutable published analysis outputs with run metadata. Immutability rule (v0.1): results are strictly write-once. Once a result is published at `/results/<result_id>/...`, its contents MUST NOT change. Deprecation MUST be represented by publishing a new result (and optionally updating status in `/registry/results` via a new registry version) and/or by updating canonical pointers.

Stored under:
- `/results/<result_id>/...`

Each result MUST have:
- `result_id` present in `/registry/results`
- `result_type` (e.g., `snv_clustering`, `clone_tree`, `cn_calls`, `mixture_model_fit`)
- `inputs` (references to registry IDs + measurement names + links used)
- `run_metadata` (who/when/how; parameters; software)

---

## 7. Canonical pointers

Canonical pointers tell users “what to use” without mutating underlying data.

Stored under:
- `/canonical/<name>`

Each canonical pointer MUST be a small record containing:
- `target_kind` (string; `registry` or `result`)
- `target_id` (string)
- `target_path` (string)
- `designated_at` (string; ISO-8601 UTC)
- `designated_by` (string)

Examples:
- `/canonical/variants_primary` → points to a variant set or a specific result that defines the active variant subset
- `/canonical/segmentation_bulk_primary`
- `/canonical/clone_tree_primary`

---

## 8. Provenance

Provenance is append-only and MUST support reproducibility.

Required tables:
- `/provenance/runs` — run records (parameters, software versions, inputs)
- `/provenance/changes` — dataset/group modifications (create/append/designate/deprecate)

For v0.1, every append-by-block operation MUST record one change row per block written, including the measurement name, block_id, grow_axis, and the affected ID range.

---

## 9. Open design decisions (to resolve next)

These are the next concrete spec decisions to lock down.

1) **Registry storage format (RESOLVED for v0.1)**
   - Registries are stored as compound (row-wise) HDF5 datasets at `/registry/<entity_name>`.

2) **ID format (RESOLVED for v0.1)**
   - All entity IDs (`*_id`) are UUIDv4 strings.
   - Human-readable identifiers belong in `label` (and optional `external_id`).

3) **Immutability enforcement mechanism (RESOLVED for v0.1)**
   - Hybrid: write-once objects + append-by-block for growing tensors.
   - Existing datasets MUST NOT be modified in-place; changes are represented as new versions or new blocks.

4) **Concurrency model**
   - single-writer lock strategy and how it’s implemented (file lock + lock record).

5) **Legacy schema migration**
   - how to migrate the current v0 layout (`/SNV`, `/sample`, `/copy_numbers`) into this target layout.

---