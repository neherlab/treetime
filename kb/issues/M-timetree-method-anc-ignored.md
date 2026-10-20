# --method-anc ignored in timetree

The ancestral reconstruction method is always marginal in the timetree pipeline, regardless of the `--method-anc` value. `--method-anc=parsimony` is accepted by clap but has no effect.

## Root cause

The `method_anc` field of `TreetimeTimetreeArgsRaw` [packages/app-commands/src/commands/timetree/args.rs#L358-L360](../../packages/app-commands/src/commands/timetree/args.rs#L358-L360) is never read. The timetree pipeline builds marginal partitions in `fn initialize_branch_model()` [packages/treetime/src/timetree/pipeline.rs#L305-L330](../../packages/treetime/src/timetree/pipeline.rs#L305-L330), updates them in each refinement round through `BranchModel::marginal_update()` [packages/treetime/src/timetree/branch_model.rs#L21-L28](../../packages/treetime/src/timetree/branch_model.rs#L21-L28), and ends with `fn final_marginal_round()` [packages/treetime/src/timetree/round.rs#L81-L92](../../packages/treetime/src/timetree/round.rs#L81-L92). No step depends on the method. The app shows the setting as "Ancestral method", so users can select a value that has no effect.

The ancestral command implements method selection correctly.

A related Python report found that a marginal-sequence option was translated correctly for the first ancestral reconstruction but omitted from a later reconstruction call, causing method selection to change during the run [[issue](https://github.com/neherlab/treetime/issues/601)]. The Rust defect has a different call path and always ignores `--method-anc`; the public report is regression context rather than a shared-root-cause claim.

## Fix direction

Implement parsimony reconstruction in timetree, or remove the flag from the timetree arguments.

## Related issues

- [M-cli-flags-parsed-but-ignored.md](M-cli-flags-parsed-but-ignored.md) lists other unused flags in
  the same pipeline
