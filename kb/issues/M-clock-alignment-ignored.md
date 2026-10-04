# `clock` ignores the alignment

`treetime clock` accepts `--aln`/`--alignment` and never reads it. v1 `clock` always works on the input branch lengths of the tree: it fits the root-to-tip regression, filters outliers, and reroots on those lengths. The alignment changes nothing, and `--covariation` still requires `--sequence-length`, although the help text of `--sequence-length` says it is "Not required if alignment is provided" ([packages/app-commands/src/commands/clock/args.rs#L185-L187](../../packages/app-commands/src/commands/clock/args.rs#L185-L187)).

The flag has existed since the first v1 command skeleton, which copied the v0 argument list. No commit has wired it, and no decision drops it: it is an unported v0 feature, not an intentional change.

## What v0 does with the alignment

v0 `estimate_clock_model()` builds a `TreeTime` object with the alignment ([packages/legacy/treetime/treetime/wrappers.py#L951-L960](../../packages/legacy/treetime/treetime/wrappers.py#L951-L960)). The alignment has two uses:

- The sequence length for the covariation variance model comes from the alignment when `--sequence-length` is absent ([packages/legacy/treetime/treetime/wrappers.py#L945-L947](../../packages/legacy/treetime/treetime/wrappers.py#L945-L947))
- With `--covariation` and without `--keep-root`, v0 calls `myTree.run(root='least-squares', max_iter=0, use_covariation=...)` before rerooting ([packages/legacy/treetime/treetime/wrappers.py#L979-L981](../../packages/legacy/treetime/treetime/wrappers.py#L979-L981)). `TreeTime._run()` reconstructs ancestral sequences and, in `joint` branch-length mode, re-estimates the branch lengths from them ([packages/legacy/treetime/treetime/treetime.py#L234-L243](../../packages/legacy/treetime/treetime/treetime.py#L234-L243)). v0 selects `joint` automatically when the longest input branch is at most 0.1, and `input` otherwise ([packages/legacy/treetime/treetime/treetime.py#L427-L455](../../packages/legacy/treetime/treetime/treetime.py#L427-L455)). The covariation regression and the root position then use the sequence-based branch lengths

### Random step in the v0 path

Every reconstruction inside `TreeTime._run()` passes `sample_from_profile='root'` ([packages/legacy/treetime/treetime/treetime.py#L211-L217](../../packages/legacy/treetime/treetime/treetime.py#L211-L217)). The root sequence is drawn from its posterior profile, one uniform draw per alignment column ([packages/legacy/treetime/treetime/treeanc.py#L1023-L1033](../../packages/legacy/treetime/treetime/treeanc.py#L1023-L1033), [packages/legacy/treetime/treetime/seq_utils.py#L266-L271](../../packages/legacy/treetime/treetime/seq_utils.py#L266-L271)). The sampled root changes the mutations on the branches below the root, and in `joint` mode their branch lengths. This draw is the only random step of v0 `clock`, and the only use of its `--rng-seed`.

## Impact

- `clock --aln=<file>` gives the same result as `clock` without an alignment, with no warning
- `clock --covariation` cannot take the sequence length from the alignment
- On trees with short branches, v1 `clock` regresses on the input branch lengths where v0 regresses on branch lengths re-estimated from the sequences, so the clock rate, the root position, and the outlier set can differ from v0

## Potential solutions

- Port the v0 path: reconstruct ancestral sequences and optimize branch lengths before the regression, with the building blocks that `optimize` and `timetree` already use. The root sampling is a random step, so the port must follow [kb/decisions/cli-seed-on-commands-with-random-steps.md](../decisions/cli-seed-on-commands-with-random-steps.md): `clock` gains `--seed` through `SeedArgs`, resolves one seed per run (drawn and logged when `--seed` is absent), and passes it to the reconstruction. Taking the most likely root state instead of sampling avoids the random step, but is a v0 divergence that needs approval
- Take only the sequence length from the alignment, and leave branch-length re-estimation to `timetree`; this is a v0 divergence that needs approval
- Remove `--aln` from `clock` and record a decision that `clock` works on tree branch lengths only

## Related issues

- [M-cli-flags-parsed-but-ignored.md](M-cli-flags-parsed-but-ignored.md) lists the other ignored `clock` flags; `--model`, `--branch-length-mode`, `--method-anc`, and `--prune-short` only matter once the alignment is used
- [M-timetree-root-sequence-not-sampled.md](M-timetree-root-sequence-not-sampled.md): the same v0 root sampling in `timetree`
- [M-clock-covariation-variance-diverges-from-v0.md](M-clock-covariation-variance-diverges-from-v0.md) records the sequence-length mismatch of `--covariation`
