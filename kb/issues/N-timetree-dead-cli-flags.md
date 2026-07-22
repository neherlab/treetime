# Dead CLI flags in timetree

Nine flags are parsed by clap but never read in the timetree pipeline:

| Flag                       | Notes                                    |
| -------------------------- | ---------------------------------------- |
| `--keep-polytomies`        | Never read                               |
| `--tip-labels`             | Never read                               |
| `--no-tip-labels`          | Never read                               |
| `--n-iqd`                  | Never read                               |
| `--vcf-reference`          | Never read                               |
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

## Potential solutions

- O1. Implement the documented behavior of a flag and add an end-to-end parse/use test.
- O2. Remove a flag and every generated-help/reference mention when no supported behavior exists.

## Recommendation

Trace each flag independently against v0 behavior and its owning feature issue, then select O1 or O2 per flag. Do not bundle all nine flags into one implementation ticket: their scientific meaning, input requirements, and parity constraints are independent.

## Ticket readiness

No aggregate ticket is ready. Create one focused ticket only after the disposition of its individual flag is determined.

## Related issues

- [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md) `--method-anc` also dead
