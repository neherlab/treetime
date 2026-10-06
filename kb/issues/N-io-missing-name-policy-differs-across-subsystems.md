# Missing-name policy differs across subsystems

FASTA records, metadata dates, and mugration traits pair with tree nodes through one shared function, `fn pair_by_name()` [packages/treetime-graph/src/pair_by_name.rs#L7](../../packages/treetime-graph/src/pair_by_name.rs#L7). It applies the same rules everywhere: one name index built while reading the entries, the first entry of a repeated name wins, and the entries without a matching node are returned for reporting. What each subsystem does with a node that has no entry still differs:

| Subsystem | Call site | Node without an entry | Entry without a node |
| --- | --- | --- | --- |
| FASTA (`ancestral`, `homoplasy`) | `fn read_nwk_fasta()` [packages/app-commands/src/commands/shared/sequence_inputs.rs](../../packages/app-commands/src/commands/shared/sequence_inputs.rs) | warning per leaf, filled with unknown characters; error above one third without `--ignore-missing-alns` | kept for the alignment length and mask, no message |
| FASTA (`timetree`, `optimize`, `prune`) | `fn pair_alignment()` [packages/app-commands/src/commands/shared/alignment.rs#L47](../../packages/app-commands/src/commands/shared/alignment.rs#L47) | error at partition creation for the first leaf | no message |
| Dates (`clock`, `timetree`) | `fn read_input_dates()` [packages/app-commands/src/commands/shared/dates_input.rs](../../packages/app-commands/src/commands/shared/dates_input.rs) | counted in the date summary; error below three dated leaves | warning with up to 10 names |
| Mugration | `fn pair_traits()` [packages/app-commands/src/commands/mugration/run.rs](../../packages/app-commands/src/commands/mugration/run.rs) | error listing every missing leaf | warning with up to 10 names |

## Impact

- The same input problem gives an error in one command and a warning in another
- Mismatch reports have different forms: one line per leaf, a capped name list, or a full list in an error

## Related

- [M-io-sequence-name-matching-unreliable.md](M-io-sequence-name-matching-unreliable.md) -- user-facing name-matching defects
- [M-mugration-tree-leaves-missing-from-metadata-rejected.md](M-mugration-tree-leaves-missing-from-metadata-rejected.md) -- the strictest of these policies
- Design: [kb/proposals/input-name-matching-validation.md](../proposals/input-name-matching-validation.md) (Architecture section, Option B: shared function)
