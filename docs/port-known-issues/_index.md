# Known Issues

Bugs, limitations, scientific and numerical issues in v1.

Known issues are distinct from deviations (intentional behavioral differences
from v0) and unimplemented features (v0 capability not yet ported).

## Filename convention

Files are prefixed with a severity letter so that `ls` sorts high-severity first:

| Prefix | Severity   | Criteria                                    |
| ------ | ---------- | ------------------------------------------- |
| `C-`   | Critical   | Data corruption, security, blocks all usage |
| `H-`   | High       | Blocks correct results or causes panics     |
| `M-`   | Medium     | Wrong results under specific conditions     |
| `N-`   | Negligible | Edge cases, missing guards, weak assertions |

The letters C < H < M < N sort alphabetically in severity order. Domain prefix
(`ancestral-`, `clock-`, `timetree-`, etc.) follows the severity letter, grouping
related issues within each tier.

Issue name in the summary table must match the H1 heading in the linked file
exactly.

## Summary

| Severity   | Scope        | Issue                                                                                                                          |
| ---------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------ |
| High       | Timetree     | [Coalescent contributions use TBP coordinates, backward pass uses calendar time](H-timetree-coalescent-coordinate-mismatch.md) |
| High       | Timetree     | [Gap character not handled in alphabet](H-timetree-gap-alphabet.md)                                                            |
| High       | Timetree     | [Zero-length branches cause panic](H-timetree-zero-length-panic.md)                                                            |
| Medium     | Ancestral    | [Dense-sparse log-likelihood divergence](M-ancestral-dense-sparse-divergence.md)                                               |
| Medium     | Ancestral    | [Sparse root invariance violation](M-ancestral-sparse-root-invariance.md)                                                      |
| Medium     | Ancestral    | [Sparse variable-site alphabet mismatch](M-ancestral-sparse-alphabet-mismatch.md)                                              |
| Medium     | Ancestral    | [Marginal reconstruction uses plain probability space](M-ancestral-marginal-probability-space.md)                              |
| Medium     | Clock        | [Clock covariation overdispersion hardcoded](M-clock-covariation-overdispersion.md)                                            |
| Medium     | Timetree     | [Date column header matching breaks on hash](M-timetree-date-header-hash.md)                                                   |
| Medium     | Timetree     | [GTR model selection not implemented](M-timetree-gtr-selection.md)                                                             |
| Medium     | Timetree     | [Polytomy zero-branch penalty differs from v0](M-timetree-polytomy-zero-branch-penalty.md)                                     |
| Medium     | Timetree     | [Positional likelihood metric differs from v0](M-timetree-positional-likelihood-metric.md)                                     |
| Medium     | Timetree     | [Skyline coalescent uses Nelder-Mead instead of SLSQP](M-timetree-skyline-nelder-mead-optimizer.md)                            |
| Negligible | Ancestral    | [Dense backward pass produces NaN for all-zero probability rows](N-ancestral-dense-normalize-log-nan.md)                       |
| Negligible | Clock        | [assign_dates fails for small trees](N-clock-small-trees.md)                                                                   |
| Negligible | Core         | [Zero branch length clamping](N-core-branch-length-clamping.md)                                                                |
| Negligible | Distribution | [Formula discretization errors silently swallowed](N-distribution-formula-silent-discretization.md)                            |
