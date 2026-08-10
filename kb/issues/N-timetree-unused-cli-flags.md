# Unused CLI flags in timetree

Eight flags are parsed by clap but never read in the timetree pipeline:

| Flag                       | Notes                                    |
| -------------------------- | ---------------------------------------- |
| `--keep-polytomies`        | Never read                               |
| `--tip-labels`             | Never read                               |
| `--no-tip-labels`          | Never read                               |
| `--n-iqd`                  | Never read                               |
| `--vcf-reference`          | Never read                               |
| `--reconstruct-tip-states` | Never read                               |
| `--report-ambiguous`       | Never read                               |
| `--model-params`           | Never read (renamed from `--gtr-params`) |

## Potential solutions

- O1. Implement the documented behavior of a flag and add an end-to-end parse/use test.
- O2. Remove a flag and every generated-help/reference mention when no supported behavior exists.

## Recommendation

Trace each flag independently against v0 behavior and its owning feature issue, then select O1 or O2 per flag. Do not bundle all eight flags into one implementation ticket: their scientific meaning, input requirements, and parity constraints are independent.

## Ticket readiness

No aggregate ticket is ready. Create one focused ticket only after the disposition of its individual flag is determined.

## Related issues

- [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md) tracks another unused argument
