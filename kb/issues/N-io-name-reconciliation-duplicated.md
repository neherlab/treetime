# Name reconciliation logic duplicated across subsystems

FASTA, dates, and mugration each reimplement tree-to-external-data name reconciliation (collect names, diff, report mismatches) at different quality levels:

| Subsystem            | Call site                      | Lookup                      | Batch errors        | Dedup |
| -------------------- | ------------------------------ | --------------------------- | ------------------- | ----- |
| FASTA (fitch)        | `fitch.rs:54-101`              | linear scan O(n\*m)         | first miss aborts   | no    |
| FASTA (marginal)     | `marginal_dense.rs:319-351`    | linear scan O(n\*m)         | first miss aborts   | no    |
| FASTA (attach)       | `attach.rs:28-89`              | `BTreeSet` O(n log n)       | warnings + count    | no    |
| Dates                | `date_constraints.rs:21-90`    | `BTreeMap::get` O(log n)    | warnings            | n/a   |
| Mugration (validate) | `marginal_discrete.rs:214-247` | `IndexSet::difference` O(n) | bidirectional batch | n/a   |
| Mugration (attach)   | `marginal_discrete.rs:47-90`   | `BTreeMap::get` O(log n)    | via validate above  | n/a   |

The mugration validator is the best implementation (bidirectional batch via `IndexSet::difference`). The FASTA sites are the worst (linear scan, first-miss abort, no dedup). Bug fixes and improvements must be applied to each site independently.

## Impact

- Inconsistent error quality across commands (mugration reports all mismatches; ancestral reports one and aborts)
- Duplicate detection absent from FASTA sites (silent wrong-data attachment)
- Future improvements (fuzzy suggestions, normalization) must be applied per site or left incomplete

## Related

- [M-io-sequence-name-matching-unreliable.md](M-io-sequence-name-matching-unreliable.md) -- user-facing defects caused by the poor reimplementations
- [M-io-sequence-attachment-quadratic.md](M-io-sequence-attachment-quadratic.md) -- O(n\*m) performance in the FASTA sites
- Design: [kb/proposals/input-name-matching-validation.md](../proposals/input-name-matching-validation.md) (Architecture section, Option B: shared function)
