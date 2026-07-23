# Sequence partitions test-support module is misspelled

[`test_scripts/lib/sequence_partions.py`](../../test_scripts/lib/sequence_partions.py) omits the second `ti` in `partitions`, and the package export preserves the misspelled import path [`test_scripts/lib/__init__.py#L9`](../../test_scripts/lib/__init__.py#L9).

## Impact

Searches for `partition` miss the module, and imports preserve an accidental spelling as if it were project vocabulary. The file is test support, so the defect affects maintenance and discoverability rather than runtime behavior.

## Required change

Rename the file to `sequence_partitions.py` with `mv`, update its package import and all consumers, and retain `SeqPartition` unless a separate naming review changes the type.

## Validation

- Python import checks cover direct and package-export imports.
- Fitch and dense ancestral test scripts still resolve `SeqPartition`.
- Repository search finds no `sequence_partions` references.
