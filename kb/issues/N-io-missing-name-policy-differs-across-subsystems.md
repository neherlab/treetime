# Missing-name policy differs across subsystems

FASTA records, metadata dates, and mugration traits pair with tree nodes through one shared function, `fn pair_by_name()` [packages/treetime-graph/src/pair_by_name.rs#L8](../../packages/treetime-graph/src/pair_by_name.rs#L8). It applies the same rules everywhere: one name index built while reading the entries, the first entry of a repeated name wins, and the entries without a matching node are returned for reporting. What each subsystem does with a node that has no entry still differs:

| Subsystem | Call site | Node without an entry | Entry without a node |
| --- | --- | --- | --- |
| FASTA (`ancestral`, `homoplasy`, `optimize`, `prune`, `timetree`) | `fn read_sequence_source()` [packages/app-commands/src/commands/shared/sequence_inputs.rs#L84-L124](../../packages/app-commands/src/commands/shared/sequence_inputs.rs#L84-L124), which calls `fn pair_alignment()` [packages/app-commands/src/commands/shared/alignment.rs#L73](../../packages/app-commands/src/commands/shared/alignment.rs#L73) | warning per leaf, filled with unknown characters; error above one third, which only `ancestral` and `homoplasy` can override with `--ignore-missing-alns` | kept for the alignment length check, left out of the mask, no message |
| Dates (`clock`, `timetree`) | `fn read_input_dates()` [packages/app-commands/src/commands/shared/dates_input.rs#L17](../../packages/app-commands/src/commands/shared/dates_input.rs#L17) | counted in the date summary; error below three dated leaves | warning with up to 10 names |
| Mugration | `fn pair_traits()` [packages/app-commands/src/commands/mugration/run.rs#L143](../../packages/app-commands/src/commands/mugration/run.rs#L143), which calls `fn pair_trait_rows()` [packages/treetime/src/mugration/trait_pairing.rs#L7](../../packages/treetime/src/mugration/trait_pairing.rs#L7) | warning per leaf and a count warning (`fn warn_leaves_without_data()`, shared with FASTA), filled with the unknown observation; no limit, as in v0 | warning with up to 10 names |

## Impact

- The same input problem gives an error in one command and a warning in another
- Mismatch reports have different forms: one line per leaf, a capped name list, or a full list in an error

> [!IMPORTANT]
> **Decision required.** The target policy is not approved. A missing leaf has a different meaning per subsystem: a leaf without a sequence or a trait value becomes missing data (only the sequence commands stop above one third), and a leaf without a date is left unconstrained. Options:
>
> - One policy for every subsystem, including one report form and one override flag
> - Subsystem-specific handling of missing nodes, with one shared report form for mismatched names

## Related

- [N-io-sequence-name-matching-unreliable.md](N-io-sequence-name-matching-unreliable.md) -- user-facing name-matching defects
- Design: [kb/proposals/input-name-matching-validation.md](../proposals/input-name-matching-validation.md) (current state of the pairing, suggestions, normalization)
