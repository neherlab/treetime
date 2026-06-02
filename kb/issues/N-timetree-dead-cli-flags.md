# Dead CLI flags in timetree

10 flags are parsed by clap but never read in the timetree pipeline:

| Flag                       | Notes                                    |
| -------------------------- | ---------------------------------------- |
| `--keep-polytomies`        | Never read                               |
| `--tip-labels`             | Never read                               |
| `--no-tip-labels`          | Never read                               |
| `--n-iqd`                  | Never read                               |
| `--vcf-reference`          | Never read                               |
| `--zero-based`             | Never read                               |
| `--reconstruct-tip-states` | Never read                               |
| `--report-ambiguous`       | Never read                               |
| `--seed`                   | Never read                               |
| `--model-params`           | Never read (renamed from `--gtr-params`) |

Removed: `--aa` (redundant with `--alphabet`, dropped in CLI args unification). Wired: `--reroot` and `--reroot-tips` now pass an explicit root selection spec into timetree rerooting.

Flags that are already wired and not part of this issue:

- `--confidence` promotes `time_marginal` from `never` to `only-final` when paired with `--covariation` or `--clock-std-dev`
- `--covariation` drives GLS clock regression params in the timetree pipeline
- `--clock-std-dev` provides user-specified rate standard deviation for rate susceptibility
- `--tip-slack` used in covariation variance computation
- `--time-marginal=always` already triggers confidence-interval extraction, so it is intentionally excluded from the dead-flag list

## Related issues

- [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md) `--method-anc` also dead
