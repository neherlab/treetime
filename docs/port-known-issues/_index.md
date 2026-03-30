# Known Issues

All unfinished v1 work items: bugs, missing features, missing output formats, dead CLI flags, stubs, unimplemented algorithms, and unintentional behavioral differences from v0.

Distinct from:

- [intentional changes](../port-intentional-changes/_index.md) - deliberate deviations from v0
- [proposals](../port-proposals/_index.md) - newly proposed features, not v0 ports

## Scope

Known issues track everything that prevents v1 from being a complete, correct implementation. This includes:

- **Bugs** - crashes, wrong results, numerical instability
- **Missing v0 features** - v0 capabilities not yet ported (algorithms, CLI flags, output formats)
- **Missing design features** - capabilities specified in `docs/algorithms/` design documents but not yet implemented in v0 or v1 (e.g., per-site rate variation)
- **Dead CLI flags** - flags parsed by clap but not wired to any runtime behavior
- **Stubs** - `todo!()`, `unimplemented!()`, functions that return placeholder values
- **Behavioral differences** - unintentional divergences from v0 results (distinct from `port-intentional-changes/`)

When a feature appears in `docs/algorithms/` design documents, it is a work item even if v0 does not implement it. Design documents are specifications, not just v0 descriptions.

## Test matrix

[`_test-matrix.md`](_test-matrix.md) tracks systematic v0/v1 comparison testing:

- Rows: exact CLI args (e.g., `--coalescent=10 --time-marginal`)
- Columns: datasets
- Cells: test results (`OK`, `CRASH-<id>`, `HANG`, `-` for untested)

Issues discovered during testing are documented here and linked from the matrix.

## Filename convention

Files are prefixed with a severity letter so that `ls` sorts high-severity first:

| Prefix | Severity   | Criteria                                                                        |
| ------ | ---------- | ------------------------------------------------------------------------------- |
| `C-`   | Critical   | Data corruption, security, blocks all usage                                     |
| `H-`   | High       | Crashes, panics, blocks correct results, or missing standard expected feature   |
| `M-`   | Medium     | Wrong results under specific conditions, or missing feature affecting workflows |
| `N-`   | Negligible | Edge cases, niche missing features, weak assertions, cosmetic gaps              |

The letters C < H < M < N sort alphabetically in severity order. Domain prefix
(`ancestral-`, `clock-`, `timetree-`, etc.) follows the severity letter, grouping
related issues within each tier.

Severity for missing features: a missing feature that most users expect (e.g., a standard phylogenetic capability) is High. A missing feature that affects specific workflows is Medium. A niche or rarely used missing feature is Negligible.

Issue name in the summary table must match the H1 heading in the linked file
exactly.

## Summary

| Severity   | Scope        | Issue                                                                                                                           |
| ---------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------- |
| ~~High~~   | ~~Timetree~~ | ~~Timetree crashes on zero-length branches with grid spacing error~~ **FIXED**                                                  |
| Medium     | Timetree     | [Golden master runner missing internal node times](M-timetree-gm-runner-missing-internal-times.md)                              |
| ~~High~~   | ~~Clock~~    | ~~Clock fails with negative rate before filtering outliers~~ **FIXED**                                                          |
| Medium     | Clock        | [Clock filter residual parity](M-clock-filter-residual-parity.md)                                                               |
| Negligible | Clock        | [Clock regression all-negative-rate divergence](N-clock-regression-all-negative-rate.md)                                        |
| Negligible | Clock        | [Clock chisq near-zero determinant](N-clock-chisq-near-zero-determinant.md)                                                     |
| Medium     | Ancestral    | [Dense-sparse log-likelihood divergence](M-ancestral-dense-sparse-divergence.md)                                                |
| Medium     | Ancestral    | [Marginal reconstruction uses plain probability space](M-ancestral-marginal-probability-space.md)                               |
| Medium     | Ancestral    | [Sparse root invariance violation](M-ancestral-sparse-root-invariance.md)                                                       |
| Medium     | Ancestral    | [Sparse variable-site alphabet mismatch](M-ancestral-sparse-alphabet-mismatch.md)                                               |
| Medium     | Core         | [Branch mutations have no unified API across partition types](M-core-branch-mutations-no-unified-api.md)                        |
| Medium     | Core         | [Dummy GTR initialization pattern across commands](M-core-dummy-gtr-initialization.md)                                          |
| Medium     | Clock        | [Clock covariation overdispersion hardcoded](M-clock-covariation-overdispersion.md)                                             |
| Medium     | Dates        | [Column auto-detection gaps in CSV readers](M-dates-column-auto-detection-gaps.md)                                              |
| Medium     | GTR          | [Per-site rate variation not implemented](M-gtr-per-site-rate-variation.md)                                                     |
| Medium     | Mugration    | [Iterative GTR inference not implemented for mugration](M-mugration-iterative-gtr.md)                                           |
| Low        | GTR          | [Site-specific GTR partition integration pending](L-gtr-site-specific-partition-integration.md)                                 |
| Low        | GTR          | [Site-specific GTR end-to-end inference test pending](L-gtr-site-specific-e2e-inference-test.md)                                |
| Medium     | Optimize     | [Optimize loop does not prune zero-length branches or resolve polytomies](M-optimize-no-topology-cleanup-in-loop.md)            |
| ~~Medium~~ | ~~Prune~~    | ~~Edge collapse uses union instead of composition for mutations~~ **FIXED**                                                     |
| Medium     | Prune        | [Prune command applies merge before prune, should be reversed](M-prune-wrong-operation-order.md)                                |
| Low        | Optimize     | [Optimize loop not extracted from I/O wrapper](L-optimize-loop-not-extracted.md)                                                |
| Medium     | Timetree     | [Timetree skips initial ML branch length optimization before time inference](M-timetree-missing-initial-branch-optimization.md) |
| Medium     | Timetree     | [--aln flag silently ignored](M-timetree-aln-flag-ignored.md)                                                                   |
| Medium     | Timetree     | [Coalescent backward pass missing leaf and root contributions](M-timetree-coalescent-missing-leaf-and-root-contributions.md)    |
| Medium     | Timetree     | [Coalescent and skyline CI excludes internal nodes](M-timetree-coalescent-ci-excludes-internal.md)                              |
| ~~Medium~~ | ~~Timetree~~ | ~~Coalescent likelihood always returns None~~ **FIXED**                                                                         |
| Medium     | Timetree     | [--coalescent-opt alone skips initial Tc pass](M-timetree-coalescent-opt-skips-initial.md)                                      |
| Medium     | Timetree     | [Date column header matching breaks on hash](M-timetree-date-header-hash.md)                                                    |
| Medium     | Timetree     | [Internal node dates missing in nexus for input branch length mode](M-timetree-internal-dates-missing-input-bl.md)              |
| ~~Medium~~ | ~~Timetree~~ | ~~Internal node dates missing at scale~~ **FIXED**                                                                              |
| ~~Medium~~ | ~~Timetree~~ | ~~Internal node dates missing with bad fixed clock rate~~ **FIXED**                                                             |
| Medium     | Timetree     | [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md)                                                            |
| Medium     | Timetree     | [Nexus output missing mutation annotations](M-timetree-nexus-missing-mutations.md)                                              |
| Medium     | Timetree     | [Positional likelihood metric differs from v0](M-timetree-positional-likelihood-metric.md)                                      |
| Medium     | Timetree     | [Skyline coalescent uses Nelder-Mead instead of SLSQP](M-timetree-skyline-nelder-mead-optimizer.md)                             |
| ~~Medium~~ | ~~Timetree~~ | ~~--time-marginal=always has no effect~~ **FIXED**                                                                              |
| Medium     | Timetree     | [Marginal dense golden master node key mismatch on ebola_20](M-timetree-dense-golden-master-node-mismatch.md)                   |
| Medium     | Timetree     | [Marginal dense timetree inference disproportionately slow for mpox dataset](M-timetree-marginal-dense-mpox-slow.md)            |
| Negligible | Ancestral    | [Sparse marginal passes still use remove/insert pattern](N-ancestral-sparse-remove-insert-pattern.md)                           |
| Negligible | Ancestral    | [Dense normalize_inplace produces NaN for all-zero probability rows](N-dense-normalize-inplace-zero-row.md)                     |
| Negligible | Clock        | [assign_dates fails for small trees](N-clock-small-trees.md)                                                                    |
| Negligible | Core         | [Zero branch length clamping](N-core-branch-length-clamping.md)                                                                 |
| Negligible | Dates        | [Imprecise date upper bound not capped at present](N-dates-imprecise-upper-bound-not-capped.md)                                 |
| Negligible | Distribution | [Formula discretization errors silently swallowed](N-distribution-formula-silent-discretization.md)                             |
| Negligible | I/O          | [Multi-segment genome input not wired](N-io-multi-segment-genome-input.md)                                                      |
| Negligible | Optimize     | [initial_guess_mixed allocates Vec\<Sub\> per edge for count only](L-optimize-initial-guess-alloc.md)                           |
| Negligible | Optimize     | [Dense/sparse equivalence test bounds undocumented](N-optimize-equivalence-bounds-undocumented.md)                              |
| Negligible | Optimize     | [Optimizer evaluation functions omit mu scaling factor](N-optimize-mu-scaling-omitted-in-evaluation.md)                         |
| Negligible | Optimize     | [Optimize command accepts only a single alignment](N-optimize-multi-alignment-input.md)                                         |
| Negligible | Optimize     | [Dense/sparse equivalence validity tests silently skip None branch lengths](N-optimize-validity-test-silent-none-skip.md)       |
| Negligible | Optimize     | [Branch length likelihood does not account for indels](N-optimize-indel-contribution-to-likelihood.md)                          |
| Negligible | Optimize     | [Initial branch length guess always overwrites input values](N-optimize-initial-guess-not-optional.md)                          |
| Negligible | Timetree     | [--dates not required, misleading error when omitted](N-timetree-dates-not-required.md)                                         |
| Negligible | Timetree     | [Dead CLI flags in timetree](N-timetree-dead-cli-flags.md)                                                                      |
| Negligible | Timetree     | [gtr.json missing for --branch-length-mode=input](N-timetree-gtr-json-missing-input-bl.md)                                      |
| Negligible | Timetree     | [timetree.json missing coalescent and skyline parameters](N-timetree-json-missing-coalescent.md)                                |
| Negligible | Timetree     | [Missing output files compared to v0](N-timetree-missing-output-files.md)                                                       |
| Negligible | Timetree     | [Missing skyline output files](N-timetree-missing-skyline-output.md)                                                            |
| Negligible | Timetree     | [--n-branches-posterior panics with todo!()](N-timetree-n-branches-posterior-unimplemented.md)                                  |
| Negligible | Timetree     | [Negative coalescent Tc accepted without validation](N-timetree-negative-coalescent-tc.md)                                      |
| Negligible | Timetree     | [write_node_dates() is a todo!() stub](N-timetree-node-dates-output-unimplemented.md)                                           |
| Negligible | Timetree     | [--plot-rtt and --plot-tree typed as Option\<usize\>](N-timetree-plot-arg-type.md)                                              |
| Negligible | Timetree     | [--plot-rtt and --plot-tree panic with todo!()](N-timetree-plot-unimplemented.md)                                               |
| Negligible | Timetree     | [--keep-polytomies and --resolve-polytomies no conflicts_with declaration](N-timetree-polytomy-flags-no-conflict.md)            |
| Negligible | Timetree     | [Stochastic polytomy resolution not implemented](N-timetree-stochastic-polytomy-unimplemented.md)                               |
| Negligible | Timetree     | [Auspice JSON output incomplete](N-timetree-auspice-json-incomplete.md)                                                         |
| Negligible | Timetree     | [Unnamed root after reroot](N-timetree-unnamed-root-after-reroot.md)                                                            |
| Negligible | Timetree     | [No tree inference from alignment](N-timetree-tree-inference-unimplemented.md)                                                  |
| Negligible | Timetree     | [Polytomy resolution numerical robustness](N-timetree-polytomy-numerical-robustness.md)                                         |
| Negligible | Timetree     | [Polytomy resolution test improvements](N-timetree-polytomy-test-improvements.md)                                               |
| Medium     | Timetree     | [Branch distribution grid uses uniform spacing](M-timetree-branch-grid-uniform-resolution.md)                                   |
| Negligible | Timetree     | [Branch grid extent uses base clock_rate, not effective rate](N-timetree-branch-grid-gamma-omitted.md)                          |

## Cross-references

- [Design documents](../algorithms/) - v1 design specifications (work items come from here too, not just v0)
- [Unimplemented algorithms](../port-algo-inventory/unimplemented.md) - v0 algorithm details for unported features
- [Feature inventory](../port-feature-inventory/_index.md) - v0/v1 feature parity tracking
- [Proposals](../port-proposals/_index.md) - new v1 features (not v0 ports, not design-doc features)
