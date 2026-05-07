# Clock command discards eleven CLI arguments without wiring

`run_clock` at [packages/treetime/src/commands/clock/run.rs#L22-L47](../../packages/treetime/src/commands/clock/run.rs#L22-L47) destructures the full `TreetimeClockArgs` and pattern-binds the following fields without using them in the function body or passing them to sub-functions:

1. `aln` -- alignment path
2. `vcf_reference` -- VCF reference
3. `gtr` -- GTR model name
4. `gtr_params` -- GTR parameters
5. `branch_length_mode` -- branch length interpretation
6. `method_anc` -- ancestral reconstruction method
7. `reroot` -- `RerootMode` enum (supported by v0, silently ignored)
8. `prune_short` -- short-branch pruning threshold
9. `plot_rtt` -- root-to-tip plot output
10. `seed` -- random seed
11. `tip_slack` -- only read inside the `covariation` branch, dead otherwise

These flags are exposed via clap and accepted by the CLI parser, so users set them expecting an effect. No warning is emitted when these flags are provided.

## Impact

- `--reroot {min_dev,least-squares,oldest}` is a standard v0 feature. Users specifying `--reroot oldest` get the default reroot behavior instead.
- `--prune-short` is silently ignored, leaving short branches in the tree.
- `--plot-rtt` is typed as `Option<usize>` (wrong type for a path-valued flag per the help text) and is never consumed.
- `--gtr`, `--gtr-params`, `--branch-length-mode`, `--method-anc` are all silently inert.

## Affected code

- Destructuring: [packages/treetime/src/commands/clock/run.rs#L22-L47](../../packages/treetime/src/commands/clock/run.rs#L22-L47)
- CLI definition: [packages/treetime/src/commands/clock/args.rs](../../packages/treetime/src/commands/clock/args.rs)

## Fix

For each flag, either wire it to the appropriate behavior or remove it from the CLI surface. If removal is chosen, document the omission in [issues/](../issues/) or [decisions/](../decisions/) as appropriate.

## Related issues

- Source: [M-clock-dead-cli-arguments.md](../issues/M-clock-dead-cli-arguments.md) -- delete after full resolution
- [clock-covariation-overdispersion-hardcoded.md](clock-covariation-overdispersion-hardcoded.md): `--tip-slack` default divergence from v0
- [N-timetree-dead-cli-flags.md](../issues/N-timetree-dead-cli-flags.md) -- similar dead-flag issue in the timetree command
