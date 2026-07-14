# Wire --method-anc flag in timetree pipeline

The ancestral reconstruction method is always marginal in the timetree pipeline, regardless of the `--method-anc` value. `--method-anc=parsimony` is accepted by clap but has no effect.

## Root cause

`args.method_anc` is declared in [`args.rs#L259`](../../packages/treetime/src/commands/timetree/args.rs#L259) but never read in the timetree pipeline. The pipeline always calls `initialize_marginal()` and `update_marginal()` regardless of the specified method.

The ancestral command implements method selection correctly.

A related Python report found that a marginal-sequence option was translated correctly for the first ancestral reconstruction but omitted from a later reconstruction call, causing method selection to change during the run [[issue](https://github.com/neherlab/treetime/issues/601)]. The Rust defect has a different call path and always ignores `--method-anc`; the public report is regression context rather than a shared-root-cause claim.

## Related issues

- Source: [kb/issues/M-timetree-method-anc-ignored.md](../issues/M-timetree-method-anc-ignored.md) -- delete after full resolution
- [Dead CLI flags in timetree](../issues/N-timetree-dead-cli-flags.md) lists other dead
  flags in the same pipeline
