# Dead CLI flags in timetree

12 flags are parsed by clap but never read in the timetree pipeline:

| Flag                       | Notes                                  |
| -------------------------- | -------------------------------------- |
| `--keep-polytomies`        | Never read                             |
| `--tip-labels`             | Never read                             |
| `--no-tip-labels`          | Never read                             |
| `--n-iqd`                  | Never read                             |
| `--vcf-reference`          | Never read                             |
| `--aa`                     | Never read                             |
| `--keep-overhangs`         | Never read                             |
| `--zero-based`             | Never read                             |
| `--reconstruct-tip-states` | Never read                             |
| `--report-ambiguous`       | Never read                             |
| `--seed`                   | Never read                             |
| `--gtr-params`             | Never read                             |
| `--reroot`                 | Accepted but always uses least-squares |

Previously dead flags now wired:

- `--confidence` promotes `time_marginal` from `never` to `only-final` when paired with `--covariation` or `--clock-std-dev`
- `--covariation` drives GLS clock regression params in the timetree pipeline
- `--clock-std-dev` provides user-specified rate standard deviation for rate susceptibility
- `--tip-slack` used in covariation variance computation

`--reroot` is accepted but `reroot_tree()` at
[`packages/treetime/src/commands/timetree/optimization/reroot.rs#L30`](../../packages/treetime/src/commands/timetree/optimization/reroot.rs#L30)
always uses `RerootParams::default()` (least-squares). The `--reroot=oldest`,
`--reroot=min-dev` modes are never dispatched.

## Related issues

- [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md) `--method-anc` also dead
- [--aln flag silently ignored](M-timetree-aln-flag-ignored.md) `--aln` also dead, causes confusing error
- [--time-marginal=always has no effect](M-timetree-time-marginal-always-ignored.md) variant exists but unchecked
