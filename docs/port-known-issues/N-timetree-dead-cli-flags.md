# Dead CLI flags in timetree

16 flags are parsed by clap but never read in the timetree pipeline:

| Flag                       | Notes                                  |
| -------------------------- | -------------------------------------- |
| `--tip-slack`              | Only functional in `clock` command     |
| `--covariation`            | Only functional in `clock` command     |
| `--keep-polytomies`        | Never read                             |
| `--tip-labels`             | Never read                             |
| `--no-tip-labels`          | Never read                             |
| `--n-iqd`                  | Never read                             |
| `--clock-std-dev`          | Never read                             |
| `--vcf-reference`          | Never read                             |
| `--aa`                     | Never read                             |
| `--keep-overhangs`         | Never read                             |
| `--zero-based`             | Never read                             |
| `--reconstruct-tip-states` | Never read                             |
| `--report-ambiguous`       | Never read                             |
| `--seed`                   | Never read                             |
| `--gtr-params`             | Never read                             |
| `--reroot`                 | Accepted but always uses least-squares |

`--tip-slack` and `--covariation` are functional in the `clock` command but not
wired into the timetree pipeline's clock regression. In v0, both affect clock
regression during timetree iterations.

`--reroot` is accepted but `reroot_tree()` at
[`reroot.rs#L30`](../../packages/treetime/src/commands/timetree/optimization/reroot.rs#L30)
always uses `RerootParams::default()` (least-squares). The `--reroot=oldest`,
`--reroot=min-dev` modes are never dispatched.

## Related issues

- [--confidence flag ignored](M-timetree-confidence-flag-ignored.md) `--confidence` also dead
- [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md) `--method-anc` also dead
- [--aln flag silently ignored](M-timetree-aln-flag-ignored.md) `--aln` also dead, causes confusing error
- [--time-marginal=always has no effect](M-timetree-time-marginal-always-ignored.md) variant exists but unchecked
