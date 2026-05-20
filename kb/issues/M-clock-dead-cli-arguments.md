# Clock command discards eight CLI arguments without wiring

`run_clock` at [packages/treetime/src/commands/clock/run.rs](../../packages/treetime/src/commands/clock/run.rs) destructures the full `TreetimeClockArgs` and pattern-binds the following fields without using them in the function body or passing them to sub-functions:

1. `aln` -- alignment path
2. `vcf_reference` -- VCF reference
3. `gtr` -- GTR model name
4. `gtr_params` -- GTR parameters
5. `branch_length_mode` -- branch length interpretation
6. `method_anc` -- ancestral reconstruction method
7. `prune_short` -- short-branch pruning threshold
8. `seed` -- random seed

Previously dead flags now wired or removed: `reroot` (now used in `RerootParams`), `tip_slack` (now used in variance calculation), `plot_rtt` (removed from CLI).

These flags are exposed via clap and accepted by the CLI parser, so users set them expecting an effect. No warning is emitted when these flags are provided.

## Impact

- `--prune-short` is silently ignored, leaving short branches in the tree.
- `--gtr`, `--gtr-params`, `--branch-length-mode`, `--method-anc` are all silently inert.

## Affected code

- Destructuring: [packages/treetime/src/commands/clock/run.rs#L22-L47](../../packages/treetime/src/commands/clock/run.rs#L22-L47)
- CLI definition: [packages/treetime/src/commands/clock/args.rs](../../packages/treetime/src/commands/clock/args.rs)

## Fix

For each flag, either wire it to the appropriate behavior or remove it from the CLI surface. If removal is chosen, document the omission in [issues/](.) or [decisions/](../decisions/) as appropriate.

## Related

- [M-clock-covariation-overdispersion.md](M-clock-covariation-overdispersion.md): `--tip-slack` default divergence from v0
- [N-timetree-dead-cli-flags.md](N-timetree-dead-cli-flags.md): similar dead-flag issue in the timetree command
