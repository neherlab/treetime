# Known Issues

Bugs, missing features, dead CLI flags, stubs, and unintentional differences from v0. Covers all v0 features not yet ported to v1.

Distinct from [intentional changes](../port-intentional-changes/_index.md) (deliberate deviations) and [proposals](../port-proposals/_index.md) (new v1 features not in v0).

## Test matrix

[`_test-matrix.md`](_test-matrix.md) tracks systematic v0/v1 comparison testing:

- Rows: exact CLI args (e.g., `--coalescent=10 --time-marginal`)
- Columns: datasets
- Cells: test results (`OK`, `CRASH-<id>`, `HANG`, `-` for untested)

Issues discovered during testing are documented here and linked from the matrix.

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

| Severity   | Scope        | Issue                                                                                                                  |
| ---------- | ------------ | ---------------------------------------------------------------------------------------------------------------------- |
| High       | Timetree     | [Gap character not handled in alphabet](H-timetree-gap-alphabet.md)                                                    |
| High       | Timetree     | [Marginal dense backward pass crash on ebola](H-timetree-marginal-dense-backward-crash.md)                             |
| High       | Timetree     | [--vary-rate panics with todo!()](H-timetree-vary-rate-unimplemented.md)                                               |
| High       | Clock        | [Clock fails with negative rate before filtering outliers](H-clock-negative-rate-before-filter.md)                     |
| Medium     | Ancestral    | [Dense-sparse log-likelihood divergence](M-ancestral-dense-sparse-divergence.md)                                       |
| Medium     | Ancestral    | [Marginal reconstruction uses plain probability space](M-ancestral-marginal-probability-space.md)                      |
| Medium     | Ancestral    | [Sparse root invariance violation](M-ancestral-sparse-root-invariance.md)                                              |
| Medium     | Ancestral    | [Sparse variable-site alphabet mismatch](M-ancestral-sparse-alphabet-mismatch.md)                                      |
| Medium     | Core         | [Branch mutations have no unified API across partition types](M-core-branch-mutations-no-unified-api.md)               |
| Medium     | Clock        | [Clock covariation overdispersion hardcoded](M-clock-covariation-overdispersion.md)                                    |
| Medium     | Dates        | [Column auto-detection gaps in CSV readers](M-dates-column-auto-detection-gaps.md)                                     |
| Medium     | Mugration    | [Iterative GTR inference not implemented for mugration](M-mugration-iterative-gtr.md)                                  |
| Medium     | Optimize     | [Optimize command hardcodes GTR to JC69](M-optimize-gtr-hardcoded-jc69.md)                                             |
| Medium     | Optimize     | [Branch length optimization oscillates without damping](M-optimize-oscillation-no-damping.md)                          |
| Medium     | Timetree     | [--aln flag silently ignored](M-timetree-aln-flag-ignored.md)                                                          |
| Medium     | Timetree     | [Coalescent and skyline CI excludes internal nodes](M-timetree-coalescent-ci-excludes-internal.md)                     |
| Medium     | Timetree     | [Coalescent likelihood always returns None](M-timetree-coalescent-likelihood-stub.md)                                  |
| Medium     | Timetree     | [--coalescent-opt alone skips initial Tc pass](M-timetree-coalescent-opt-skips-initial.md)                             |
| Medium     | Timetree     | [--confidence flag ignored](M-timetree-confidence-flag-ignored.md)                                                     |
| Medium     | Timetree     | [Date column header matching breaks on hash](M-timetree-date-header-hash.md)                                           |
| Medium     | Timetree     | [GTR model selection not implemented](M-timetree-gtr-selection.md)                                                     |
| Medium     | Timetree     | [Internal node dates missing in nexus for input branch length mode](M-timetree-internal-dates-missing-input-bl.md)     |
| Medium     | Timetree     | [Internal node dates missing at scale](M-timetree-internal-dates-missing-scale.md)                                     |
| Medium     | Timetree     | [Internal node dates missing with bad fixed clock rate](M-timetree-internal-dates-bad-fixed-rate.md)                   |
| Medium     | Timetree     | [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md)                                                   |
| Medium     | Timetree     | [Nexus output missing mutation annotations](M-timetree-nexus-missing-mutations.md)                                     |
| Medium     | Timetree     | [Polytomy zero-branch penalty differs from v0](M-timetree-polytomy-zero-branch-penalty.md)                             |
| Medium     | Timetree     | [Positional likelihood metric differs from v0](M-timetree-positional-likelihood-metric.md)                             |
| Medium     | Timetree     | [Skyline coalescent uses Nelder-Mead instead of SLSQP](M-timetree-skyline-nelder-mead-optimizer.md)                    |
| Medium     | Timetree     | [--time-marginal=always has no effect](M-timetree-time-marginal-always-ignored.md)                                     |
| Negligible | Ancestral    | [Dense backward pass produces NaN for all-zero probability rows](N-ancestral-dense-normalize-log-nan.md)               |
| Negligible | Clock        | [assign_dates fails for small trees](N-clock-small-trees.md)                                                           |
| Negligible | Core         | [Zero branch length clamping](N-core-branch-length-clamping.md)                                                        |
| Negligible | Dates        | [Imprecise date upper bound not capped at present](N-dates-imprecise-upper-bound-not-capped.md)                        |
| Negligible | Optimize     | [Standalone sparse optimizer skips derivative check for zero branches](N-optimize-sparse-zero-branch-no-derivative.md) |
| Negligible | Distribution | [Formula discretization errors silently swallowed](N-distribution-formula-silent-discretization.md)                    |
| Negligible | Timetree     | [--dates not required, misleading error when omitted](N-timetree-dates-not-required.md)                                |
| Negligible | Timetree     | [Dead CLI flags in timetree](N-timetree-dead-cli-flags.md)                                                             |
| Negligible | Timetree     | [gtr.json missing for --branch-length-mode=input](N-timetree-gtr-json-missing-input-bl.md)                             |
| Negligible | Timetree     | [timetree.json missing coalescent and skyline parameters](N-timetree-json-missing-coalescent.md)                       |
| Negligible | Timetree     | [Missing output files compared to v0](N-timetree-missing-output-files.md)                                              |
| Negligible | Timetree     | [Missing skyline output files](N-timetree-missing-skyline-output.md)                                                   |
| Negligible | Timetree     | [--n-branches-posterior panics with todo!()](N-timetree-n-branches-posterior-unimplemented.md)                         |
| Negligible | Timetree     | [Negative coalescent Tc accepted without validation](N-timetree-negative-coalescent-tc.md)                             |
| Negligible | Timetree     | [write_node_dates() is a todo!() stub](N-timetree-node-dates-output-unimplemented.md)                                  |
| Negligible | Timetree     | [--plot-rtt and --plot-tree typed as Option\<usize\>](N-timetree-plot-arg-type.md)                                     |
| Negligible | Timetree     | [--plot-rtt and --plot-tree panic with todo!()](N-timetree-plot-unimplemented.md)                                      |
| Negligible | Timetree     | [--keep-polytomies and --resolve-polytomies no conflicts_with declaration](N-timetree-polytomy-flags-no-conflict.md)   |
| Negligible | Timetree     | [Stochastic polytomy resolution not implemented](N-timetree-stochastic-polytomy-unimplemented.md)                      |
| Negligible | Timetree     | [No tree inference from alignment](N-timetree-tree-inference-unimplemented.md)                                         |

## Cross-references

- [Unimplemented algorithms](../port-algo-inventory/unimplemented.md) - v0 algorithm details for unported features
- [Feature inventory](../port-feature-inventory/_index.md) - v0/v1 feature parity tracking
- [Proposals](../port-proposals/_index.md) - new v1 features (not v0 ports)
