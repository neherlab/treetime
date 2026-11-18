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

## Proposed implementation plan

1. `mv test_scripts/lib/sequence_partions.py test_scripts/lib/sequence_partitions.py`.
2. In `test_scripts/lib/__init__.py#L9`, change `from .sequence_partions import *` to `from .sequence_partitions import *`. This is the only import of the module; the consumers `test_scripts/fitch.py#L5` and `test_scripts/ancestral_dense.py#L5` import `SeqPartition` through the package export `from lib import ...` and need no change.
3. Delete the stale `test_scripts/lib/__pycache__/sequence_partions.*.pyc` file if present (untracked build output).
4. No new test: the scripts have no test suite, and the rename has no behavior.

Validation:

- `rg sequence_partions -g '!kb/**'` finds nothing.
- Inside the container, `python -c 'from lib import SeqPartition'` run from `test_scripts/` succeeds, and `python -c 'import fitch'` and `python -c 'import ancestral_dense'` resolve their imports (or run both scripts on their default inputs).
