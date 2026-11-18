# Homoplasy: filled overhangs dominate the ambiguous-change lists

> [!IMPORTANT]
> **Decision required.** The ambiguous-change lists of `treetime homoplasy` count the alignment ends that gap filling turns into `N`. These columns are missing data from sequencing coverage, not changes, and they fill the top of the ranked list. Options and evidence are below.

## Symptom

`homoplasy` sorts every change involving an ambiguous character into its own lists ([kb/decisions/homoplasy-mutation-mapping-and-counting.md](../decisions/homoplasy-mutation-mapping-and-counting.md)). Before reconstruction, `--gap-fill=only-terminal` (the default, matching the v0 overhang filling) replaces the leading and trailing gaps of each sequence with `N` (`fn read_sequence_source()` [packages/app-commands/src/commands/shared/sequence_inputs.rs#L265-L267](../../packages/app-commands/src/commands/shared/sequence_inputs.rs#L265-L267)). Each filled column of a sample whose parent has a determined state is then a change to `N` on the terminal branch.

On `data/rsv/a/20`, 3028 ambiguous changes are counted. The 15 most frequent are `G1N`, `T2N`, `A3N`, and the following columns of the 5' end, each on 9 or 10 terminal branches, and the per-taxon table reports more than 1000 ambiguous changes for the samples `JF920067` and `KJ627697`:

```bash
./dev/docker/run just r treetime homoplasy --tree=data/rsv/a/20/tree.nwk --alignment=data/rsv/a/20/aln.fasta.xz --detailed -n 15 --output-all=tmp/homoplasy/rsv20
```

The substitution statistics are not affected: changes to `N` never enter them.

## Options

- O1. Keep counting every change to `N`. The lists show the coverage of each sample, and the per-taxon count of a sample with short reads is large by construction
- O2. Count only ambiguous characters that the input sequences contain, and leave out the columns that gap filling created. The pipeline would need to tell filled columns from observed `N`, for example from the gap-fill mask of each sample
- O3. Count changes to partial ambiguity codes (such as `R` or `Y`) only, and leave out every `N`. Observed `N` in the middle of a sequence then disappears from the lists as well

## Proposed implementation plan

`fn read_sequence_source()` applies `apply_gap_fill()` at the cited lines, and the homoplasy ambiguous lists live in `packages/app-commands/src/commands/homoplasy/` (`summary.rs`, `report.rs`, `result.rs`). The counts on `data/rsv/a/20` come from earlier runs and were not measured again on the current code. The plan stops at the open decision; one short plan per option:

- **O1**: delete this issue; optionally document the coverage meaning of the lists in the `--help` text of `homoplasy`.
- **O2**: make `apply_gap_fill()` return the filled ranges per sample, carry them with the alignment source to the homoplasy count, and skip changes to `N` inside a filled range. Test: a 3-tip alignment with one leading-gap sample and one observed `N`; the ambiguous list holds only the observed `N` (oracle: manual count).
- **O3**: in the homoplasy ambiguous filter, skip changes whose target is `N`. Test: same alignment with an `R` and an `N`; the list holds only the `R` change (oracle: manual count).

For O2 and O3, update [kb/decisions/homoplasy-mutation-mapping-and-counting.md](../decisions/homoplasy-mutation-mapping-and-counting.md) with approval. Validation: `just lint-rs`, `just test-rs`, and the command above, whose top 15 ambiguous changes then contain no 5' end columns.
