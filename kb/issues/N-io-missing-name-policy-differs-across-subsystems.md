# Missing-name policy differs across subsystems

FASTA records, metadata dates, and mugration traits pair with tree nodes through one shared function, `fn pair_by_name()` [packages/treetime-graph/src/pair_by_name.rs#L9](../../packages/treetime-graph/src/pair_by_name.rs#L9). It applies the same rules everywhere: one name index built while reading the entries, the first entry of a repeated name wins, and the entries without a matching node are returned for reporting. What each subsystem does with a node that has no entry still differs:

| Subsystem | Call site | Node without an entry | Entry without a node |
| --- | --- | --- | --- |
| FASTA (`ancestral`, `homoplasy`, `optimize`, `prune`, `timetree`) | `fn read_sequence_source()` [packages/app-commands/src/commands/shared/sequence_inputs.rs#L256-L296](../../packages/app-commands/src/commands/shared/sequence_inputs.rs#L256-L296), which calls `fn pair_alignment()` [packages/app-commands/src/commands/shared/alignment.rs#L60](../../packages/app-commands/src/commands/shared/alignment.rs#L60) | warning per leaf, filled with unknown characters; error above one third, which `--allow-missing-sequences` overrides | kept for the alignment length check, left out of the mask, no message |
| Dates (`clock`, `timetree`) | `fn read_input_dates()` [packages/app-commands/src/commands/shared/dates_input.rs#L16](../../packages/app-commands/src/commands/shared/dates_input.rs#L16) | counted in the date summary; error below three dated leaves | warning with up to 10 names |
| Mugration | `fn pair_traits()` [packages/app-commands/src/commands/mugration/run.rs#L149](../../packages/app-commands/src/commands/mugration/run.rs#L149), which calls `fn pair_trait_rows()` [packages/treetime/src/mugration/trait_pairing.rs#L7](../../packages/treetime/src/mugration/trait_pairing.rs#L7) | warning per leaf and a count warning (`fn warn_leaves_without_data()`, shared with FASTA), filled with the unknown observation; no limit, as in v0 | warning with up to 10 names |

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

## Proposed implementation plan

The target policy is open (Decision required above). One short plan per option; both share step 1.

1. Shared report: add one function beside `fn warn_leaves_without_data()` (`packages/treetime/src/ancestral/attach.rs`) that reports unmatched entries and nodes without data in one form (count plus up to 10 names), and call it from `fn pair_alignment()` (`packages/app-commands/src/commands/shared/alignment.rs`), `fn read_input_dates()` (`packages/app-commands/src/commands/shared/dates_input.rs`), and `fn pair_traits()` (`packages/app-commands/src/commands/mugration/run.rs`). The FASTA path then also reports entries without a node
2. Option "one policy for every subsystem": move the one-third limit and its override flag into the shared function, rename `--allow-missing-sequences` to a flag that every command with named inputs accepts, and apply it to dates and traits. Update `packages/app-cli/src/cli/v0_args.rs` only if a 0.x spelling is affected, then `just gen`
3. Option "subsystem-specific handling": keep the limit in the FASTA path and the three-dated-leaves rule in the dates path; only step 1 changes behavior
4. Tests: one case per subsystem in the `__tests__/` directory of each call site, with a tree of four leaves and inputs that miss one leaf and add one unknown name; assert the whole returned result and the error text with `assert_err!` where the policy stops the run (oracle: the approved policy text)
5. Validate: `just lint-rs`, `just test-rs`, `just generated-check` when a flag changes
