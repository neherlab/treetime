# `clock` ignores the alignment

`treetime clock` accepts `--aln`/`--alignment` and never reads it. v1 `clock` always works on the input branch lengths of the tree: it fits the root-to-tip regression, filters outliers, and reroots on those lengths. The alignment changes nothing.

The flag has existed since the first v1 command skeleton, which copied the v0 argument list. No commit has wired it, and no decision drops it: it is an unported v0 feature, not an intentional change.

## What v0 does with the alignment

v0 `estimate_clock_model()` builds a `TreeTime` object with the alignment ([packages/legacy/treetime/treetime/wrappers.py#L951-L960](../../packages/legacy/treetime/treetime/wrappers.py#L951-L960)). Both uses of the alignment belong to `--covariation`, which v1 removed ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md)):

- The sequence length for the covariation variance model comes from the alignment when `--sequence-length` is absent ([packages/legacy/treetime/treetime/wrappers.py#L945-L947](../../packages/legacy/treetime/treetime/wrappers.py#L945-L947))
- With `--covariation` and without `--keep-root`, v0 calls `myTree.run(root='least-squares', max_iter=0, use_covariation=...)` before rerooting ([packages/legacy/treetime/treetime/wrappers.py#L979-L981](../../packages/legacy/treetime/treetime/wrappers.py#L979-L981)). `TreeTime._run()` reconstructs ancestral sequences and, in `joint` branch-length mode, re-estimates the branch lengths from them ([packages/legacy/treetime/treetime/treetime.py#L234-L243](../../packages/legacy/treetime/treetime/treetime.py#L234-L243)). v0 selects `joint` automatically when the longest input branch is at most 0.1, and `input` otherwise ([packages/legacy/treetime/treetime/treetime.py#L427-L455](../../packages/legacy/treetime/treetime/treetime.py#L427-L455)). The covariation regression and the root position then use the sequence-based branch lengths

### Random step in the v0 path

Every reconstruction inside `TreeTime._run()` passes `sample_from_profile='root'` ([packages/legacy/treetime/treetime/treetime.py#L211-L217](../../packages/legacy/treetime/treetime/treetime.py#L211-L217)). The root sequence is drawn from its posterior profile, one uniform draw per alignment column ([packages/legacy/treetime/treetime/treeanc.py#L1023-L1033](../../packages/legacy/treetime/treetime/treeanc.py#L1023-L1033), [packages/legacy/treetime/treetime/seq_utils.py#L266-L271](../../packages/legacy/treetime/treetime/seq_utils.py#L266-L271)). The sampled root changes the mutations on the branches below the root, and in `joint` mode their branch lengths. This draw is the only random step of v0 `clock`, and the only use of its `--rng-seed`.

## Impact

- `clock --aln=<file>` gives the same result as `clock` without an alignment, with no warning. Only the `--help` text of `--alignment` states that `clock` does not read the alignment (`struct AlignmentArgs` [packages/app-commands/src/commands/shared/alignment.rs#L33-L37](../../packages/app-commands/src/commands/shared/alignment.rs#L33-L37)); `clock` still accepts the flag through `TreetimeClockArgsRaw` [packages/app-commands/src/commands/clock/args.rs#L70-L74](../../packages/app-commands/src/commands/clock/args.rs#L70-L74)
- Without `--covariation`, v0 does not use the alignment either; with it removed from v1, the v0 behavior that the flag served has no v1 counterpart

## Potential solutions

> [!IMPORTANT]
> **Decision required.** `clock` accepts an alignment that it never reads, and no decision says whether `clock` works on tree branch lengths only. Options:
>
> - **Remove the flag**: remove `--aln` from `clock` and record a decision that `clock` works on tree branch lengths only
> - **Use the alignment**: re-estimate the branch lengths from the sequences before the regression, as v0 does under `--covariation`, with the building blocks that `optimize` and `timetree` already use; the root sampling is a random step, so this must follow [kb/decisions/cli-seed-on-commands-with-random-steps.md](../decisions/cli-seed-on-commands-with-random-steps.md)
>
> Evidence: v0 uses the alignment only under `--covariation`, which v1 removed (see above); the second option also gives `clock` a sequence length, which [M-clock-command-reports-no-rate-uncertainty.md](M-clock-command-reports-no-rate-uncertainty.md) needs for a rate uncertainty

## Related issues

- [M-cli-flags-parsed-but-ignored.md](M-cli-flags-parsed-but-ignored.md) lists the other ignored `clock` flags; `--model`, `--branch-length-mode`, `--method-anc`, and `--prune-short` only matter once the alignment is used
- [M-timetree-root-sequence-not-sampled.md](M-timetree-root-sequence-not-sampled.md): the same v0 root sampling in `timetree`
