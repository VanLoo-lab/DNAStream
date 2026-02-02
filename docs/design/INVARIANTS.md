# DNAStream Invariants

This document defines invariants for DNAStream.
An *invariant* is a rule that must always hold for any valid DNAStream object.
Any file or operation that violates these invariants is considered invalid.

These invariants apply across all schema versions unless explicitly superseded.

⸻

##  Identity and Registries

1. Stable Entity Identity

All Entities (Sample, Variant, Segmentation, Segment, Bin, Result) SHALL have a stable, immutable identifier (*_id) that is unique within its registry.

2. Registry Immutability

Entity registries are append-only.
Entities SHALL NOT be deleted or modified in place.
Deprecation or exclusion MUST be represented via explicit metadata e.g., active=false,deprecated_at, superseded_by.

3. Identifier Reuse Prohibited

Entity identifiers SHALL NOT be reused, even if an entity is deprecated or inactive.

4. Explicit Entity Membership

Every Segment SHALL belong to exactly one Segmentation.
Every Bin SHALL belong to exactly one Segmentation.
Segments and Bins SHALL NOT exist outside a Segmentation context.

⸻

## Measurements

I5. Measurement Keying

Measurements SHALL be stored as typed tensors indexed by one or more Entity registries (e.g., Variant × Sample, Segment × Sample).

6. Axis Declaration

Every Measurement SHALL explicitly declare its indexing entities and axis order.
Measurement shapes MUST exactly match the cardinality of the declared entity axes.

7. Measurement Non-Identity

Measurements SHALL NOT have independent identity, versioning, or canonical designation.
Measurements derive meaning solely from the entities they index.

8. Measurement Immutability

Measurements representing processed observational data (e.g., read counts, BAFs, logR) SHALL be immutable after attachment.
Modified or filtered measurements MUST be represented as new measurements or views.

⸻

## Links and Referential Integrity

I9. Explicit Links

Relationships between Entities (e.g., Variant → Segment, Bin → Segment) SHALL be stored explicitly as link objects.

I10. Referential Integrity

All links SHALL reference existing entity identifiers.
Links referencing non-existent or deprecated entities SHALL be invalid.

I11. Identity-Preserving Subsetting

Subsetting operations SHALL preserve entity identifiers and maintain referential integrity across all associated Measurements, Links, and Results.

⸻

## Results

I12. Result Identity

Every Result SHALL have a stable, immutable result_id and be stored in a Result registry.

I13. Result Immutability

Results SHALL be immutable after publication.
Any modification to a Result (including parameters, structure, or interpretation) MUST create a new Result with a new result_id.

I14. Result Scope

Results SHALL reference existing Entities, Measurements, and/or Links.
Results SHALL NOT duplicate raw Measurement data unless explicitly documented and justified.

I15. Result Provenance

Every Result SHALL be associated with provenance metadata, including at minimum:
	•	run_id
	•	creation timestamp
	•	method or tool identifier
	•	relevant parameters

⸻

## Versioning and Canonical References

I16. Schema Version Declaration

Every DNAStream object SHALL declare a schema_version.
Loaders SHALL refuse to open objects with incompatible major schema versions.

I17. Canonical Reference Semantics

Canonical references SHALL be explicit pointers to immutable Entities or Results.
Changing a canonical reference SHALL NOT mutate the referenced object.

I18. Canonical Reference Uniqueness

At most one canonical reference MAY exist per reference category (e.g., canonical SNV set, canonical CN Result) within a given scope.

I19. Canonical Reference Independence

Canonical references SHALL be independent of creation time or insertion order.
“Latest” or “most recent” SHALL NOT be used as implicit canonical behavior.

⸻

## Provenance and Access Control

I20. Run Attribution

Any operation that creates or publishes Entities, Measurements, Links, or Results SHALL be associated with a run_id.

I21. Provenance Immutability

Provenance records SHALL be append-only.
Historical provenance entries SHALL NOT be modified or deleted.

I22. Controlled Write Access

DNAStream SHALL enforce multi-user safe access, supporting concurrent readers and controlled write access (e.g., single-writer semantics) to prevent data corruption.

⸻

## Validation and Failure Semantics

I23. Validation on Write

All write or modification operations SHALL validate invariants prior to completion.
Operations that violate invariants SHALL fail.

I24. Fail Loudly and Early

Invalid shapes, mismatched indices, missing identifiers, schema violations, or referential integrity errors SHALL raise explicit errors.
Silent coercion or implicit correction is prohibited.

⸻

## Backward Compatibility and Migration

I25. Backward Compatibility

Minor schema versions SHALL remain backward compatible unless explicitly documented otherwise.

I26. Explicit Migration

Breaking schema changes SHALL require explicit migration steps.
Automatic silent migration is prohibited.

⸻

Notes

These invariants define the minimum correctness guarantees of DNAStream.
Implementations MAY impose additional constraints, but SHALL NOT violate these invariants.