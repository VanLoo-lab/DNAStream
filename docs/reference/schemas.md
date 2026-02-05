# Built-in schemas

These tables are generated from the current `dnastream` code.

## Registries

### `sample`

| name | dtype | required | validator |
| --- | --- | --- | --- |
| id | object | True | None |
| label | object | True | None |
| idx | <class 'numpy.int64'> | True | None |
| active | <class 'numpy.bool'> | True | None |
| created_at | object | True | None |
| created_by | object | True | None |
| modified_at | object | True | None |
| modified_by | object | True | None |
| sample_name | object | True | str_validator |
| tissue_type | object | False | str_validator |
| organism | object | False | str_validator |
| library_strategy | object | False | str_validator |
| library_source | object | False | str_validator |
| library_selection | object | False | str_validator |
| library_layout | S10 | False | None |
| read_length | i4 | False | None |
| platform | object | False | str_validator |
| model | object | False | str_validator |
| center_name | object | False | str_validator |
| run | object | False | str_validator |
| study | object | False | str_validator |
| coverage | f4 | False | None |
| modality | object | False | str_validator |
| location | object | False | str_validator |
| bam_file_path | object | False | str_validator |
| batch_id | object | False | str_validator |
| reference_build | object | False | str_validator |
| date_of_sequencing | object | False | str_validator |

### `variant`

| name | dtype | required | validator |
| --- | --- | --- | --- |
| id | object | True | None |
| label | object | True | None |
| idx | <class 'numpy.int64'> | True | None |
| active | <class 'numpy.bool'> | True | None |
| created_at | object | True | None |
| created_by | object | True | None |
| modified_at | object | True | None |
| modified_by | object | True | None |
| chrom | object | False | str_validator |
| start_pos | i8 | False | None |
| end_pos | i8 | False | None |
| ref_allele | object | False | str_validator |
| alt_allele | object | False | str_validator |
| hugo | object | False | str_validator |
| entrez_gene_id | object | False | str_validator |
| variant_classification | object | False | str_validator |
| variant_type | object | False | str_validator |
| dbsnp_id | object | False | str_validator |
| filter | object | False | str_validator |
| info | object | False | str_validator |
| source | object | False | str_validator |

### `snp`

| name | dtype | required | validator |
| --- | --- | --- | --- |
| id | object | True | None |
| label | object | True | None |
| idx | <class 'numpy.int64'> | True | None |
| active | <class 'numpy.bool'> | True | None |
| created_at | object | True | None |
| created_by | object | True | None |
| modified_at | object | True | None |
| modified_by | object | True | None |
| chrom | S10 | True | None |
| start_pos | i8 | True | None |
| ref_allele | S10 | True | None |
| alt_allele | S10 | True | None |
| dbsnp_id | S20 | False | None |
| strand | S1 | False | None |

## Provenance

### `provenance/log`

| name | dtype | required | validator |
| --- | --- | --- | --- |
| id | object | True | None |
| timestamp | object | True | None |
| user | object | True | None |
| scope | object | True | None |
| event | object | True | validate_event |
| dataset | object | True | None |
| source | object | False | None |
| info | object | False | None |
