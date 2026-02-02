
# Registry Contract

## Location
- All registries are stored under the HDF5 group **`/registry`** to ensure they are centrally located and discoverable.
- Each registry dataset lives at the path: `/registry/<name>`


---

## Schema Requirements
- Registry creation **must** be provided a schema.
- The schema **must** define the following fields in the dataset dtype:
- `id`
- `label`
- `idx`
- `active`
- `created_at`
- `created_by`
- `modified_at`
- `modified_by`
- The schema **must** define how labels are generated and validated:
- `label_from`: tuple of field names used to construct the label
- `label_normalizer`: callable that produces a normalized label
- Registry IDs:
    - Are UUIDs
    - Are globally unique
    - Are immutable once written

---

## Mutability Rules
- Registries are **append-only**:
- Rows may never be deleted
- Rows may never be reordered
- Existing rows may only be modified via state or metadata fields (e.g., `active`). 
    - Editable metadata fields are any fields that are not in the Registry spine (see above) and not in the `label_from` tuple
    - Any metadata edits will reset the `modified_at` and `modified_by` fields


---

## State Model
Each registry entity may exist in either an active or deactive state.

---

## Label Uniqueness
- Labels **may** be duplicated across the full registry history.
- Labels **must be unique among active entities**.
- No two rows with `active=True` may share the same normalized label.

---

## Caching
- Registries may maintain in-memory caches for:
- `label ↔ id` mappings in active_only mode
- `label ↔ idx`mappings in active only mode
- `id ↔ label`
- `idx ↔ label`
- `id ↔ idx`
- `idx ↔ id`
- All caches **must be invalidated** on any operation that mutates or resizes the dataset, including:
    - `add()`
    - `activate()`
    - `deactivate()`


---

## Design Principle
**Registries provide identity, not ordering.**

- Registries manage stable identifiers and labels.
- Ordering and indexing are defined by measurement datasets, not registries.