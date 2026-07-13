# Known Issues

All open v1 work items: bugs, missing features, missing output formats, dead CLI flags, stubs, unimplemented algorithms, and unintentional behavioral differences from v0.

Distinct from:

- [decisions](../decisions/README.md) - deliberate design choices
- [proposals](../proposals/README.md) - new features pre-implementation

## Scope

Known issues track everything that prevents v1 from being a complete, correct implementation. This includes:

- Bugs: crashes, wrong results, numerical instability
- Missing features: v0 capabilities not yet implemented (algorithms, CLI flags, output formats)
- Missing design features: capabilities specified in `_raw/` design documents but not yet implemented (e.g., per-site rate variation)
- Dead CLI flags: flags parsed by clap but not wired to any runtime behavior
- Stubs: `todo!()`, `unimplemented!()`, functions that return placeholder values
- Behavioral differences: unintentional divergences from v0 results (distinct from `decisions/`)

When a feature appears in `_raw/` design documents, it is a work item even if v0 does not implement it. Design documents are specifications.

## Filename convention

Files are prefixed with a severity letter so that `ls` sorts high-severity first:

| Prefix | Severity   | Criteria                                                                        |
| ------ | ---------- | ------------------------------------------------------------------------------- |
| `H-`   | High       | Crashes, panics, blocks correct results, or missing standard expected feature   |
| `M-`   | Medium     | Wrong results under specific conditions, or missing feature affecting workflows |
| `N-`   | Negligible | Edge cases, niche missing features, weak assertions, cosmetic                   |

The letters H < M < N sort alphabetically in severity order. Domain prefix (`ancestral-`, `clock-`, `timetree-`, etc.) follows the severity letter, grouping related issues within each tier.

Severity for missing features: a missing feature that most users expect (e.g., a standard phylogenetic capability) is High. A missing feature that affects specific workflows is Medium. A niche or rarely used missing feature is Negligible.

Issue name in the summary table must match the H1 heading in the linked file exactly.

## Summary

| Severity   | Scope        | Issue                                                                                                                                               |
| ---------- | ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| High       | Clock        | [ClockSet dateless-leaf contribution biases regression](H-clock-clockset-dateless-leaf-bias.md)                                                     |
| High       | Core         | [Command modules contain shared operations that belong in domain layers](H-core-command-module-shared-ops-entanglement.md)                          |
| High       | Core         | [Library crate is not consumable by non-CLI clients](H-core-multi-client-architecture-library-purity.md)                                            |
| High       | Homoplasy    | [Homoplasy command is unimplemented](H-homoplasy-command-unimplemented.md)                                                                          |
| High       | Timetree     | [Coalescent backward pass grid explosion](H-timetree-coalescent-grid-explosion.md)                                                                  |
| High       | Timetree     | [No tree inference from alignment](H-timetree-tree-inference-unimplemented.md)                                                                      |
| Medium     | Timetree     | [Golden master runner tests missing internal node times for 5 datasets](M-timetree-gm-runner-missing-internal-times.md)                             |
| Medium     | Clock        | [Clock filter residual computation differs from v0](M-clock-filter-residual-parity.md)                                                              |
| Medium     | Clock        | [MinDev reroot uses wrong objective (EstimatedRate instead of FixedRate(0))](M-clock-mindev-wrong-objective.md)                                     |
| Negligible | Clock        | [v1 clock regression produces all-negative rates where v0 finds positive](N-clock-regression-all-negative-rate.md)                                  |
| Negligible | Clock        | [ClockSet::chisq() numerically unstable for near-zero determinant](N-clock-chisq-near-zero-determinant.md)                                          |
| Medium     | Ancestral    | [Dense-sparse log-likelihood divergence](M-ancestral-dense-sparse-divergence.md)                                                                    |
| Medium     | Ancestral    | [Marginal reconstruction uses plain probability space](M-ancestral-marginal-probability-space.md)                                                   |
| Medium     | Ancestral    | [Sparse root invariance violation](M-ancestral-sparse-root-invariance.md)                                                                           |
| Medium     | Ancestral    | [Sparse variable-site alphabet mismatch](M-ancestral-sparse-alphabet-mismatch.md)                                                                   |
| Medium     | Ancestral    | [GenBank annotation format not supported for AA reconstruction](M-ancestral-genbank-annotation-unsupported.md)                                      |
| Medium     | Clock        | [Clock covariation overdispersion hardcoded](M-clock-covariation-overdispersion.md)                                                                 |
| Medium     | Core         | [Partition creation is hardcoded per-command instead of configured](M-core-partition-init-orchestration-duplication.md)                             |
| Medium     | Dates        | [Column auto-detection gaps in CSV readers](M-dates-column-auto-detection-gaps.md)                                                                  |
| Medium     | GTR          | [Per-site rate variation not implemented](M-gtr-per-site-rate-variation.md)                                                                         |
| Medium     | I/O          | [Sequence attachment has O(n squared) complexity](M-io-sequence-attachment-quadratic.md)                                                            |
| Medium     | I/O          | [Sequence-to-node name matching is unreliable](M-io-sequence-name-matching-unreliable.md)                                                           |
| Medium     | Mugration    | [Mugration golden master parity with v0](M-mugration-iterative-gtr.md)                                                                              |
| Medium     | Optimize     | [Indel event count reduced after composition on merged edges](M-optimize-indel-event-count-after-composition.md)                                    |
| Medium     | Optimize     | [Optimizer and marginal propagation use incompatible branch-length scales when GTR mu != 1](M-optimize-gtr-mu-coordinate-mismatch.md)               |
| Negligible | GTR          | [Site-specific GTR not integrated into partition system](N-gtr-site-specific-partition-integration.md)                                              |
| Negligible | GTR          | [Site-specific GTR inference lacks end-to-end test from real tree data](N-gtr-site-specific-e2e-inference-test.md)                                  |
| Negligible | Partition    | [`infer_dense()` stub always returns false](N-representation-infer-dense-stub.md)                                                                   |
| Negligible | Partition    | [Dense and sparse partition types have structural and naming asymmetries](N-representation-dense-sparse-partition-asymmetry.md)                     |
| Medium     | Optimize     | [Bifurcating root optimization: independent per-edge vs joint arc optimization](M-optimize-root-bifurcating-independent-vs-joint.md)                |
| Medium     | Optimize     | [Optimize per-branch lengths diverge from v0 fixture beyond 10% relative tolerance](M-optimize-gm-per-branch-divergence.md)                         |
| Medium     | Optimize     | [Sparse branch optimization reads stale Fitch-era states](M-optimize-sparse-stale-fitch-states.md)                                                  |
| Medium     | Sparse       | [Sparse marginal fixed_counts zero-multiplicity produces NaN](M-sparse-marginal-zero-multiplicity-nan.md)                                           |
| Medium     | Sparse       | [Sparse marginal msg_to_child includes near-deterministic variable sites](M-sparse-marginal-quantitative-site-filtering.md)                         |
| Medium     | Sparse       | [Sparse node_reference_state ignores non-char positions](M-sparse-node-reference-state-non-char-precedence.md)                                      |
| Medium     | Distribution | [Distribution support-boundary semantics are unresolved](M-distribution-support-boundary-semantics-unresolved.md)                                  |
| Medium     | Distribution | [Mixed-support distribution operations lose exact intersection boundaries](M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md) |
| Medium     | Distribution | [Distribution product grid resolution diverges from v0](M-distribution-product-grid-resolution-diverges-from-v0.md)                                |
| Medium     | Timetree     | [Coalescent backward pass missing leaf and root contributions](M-timetree-coalescent-missing-leaf-and-root-contributions.md)                        |
| Medium     | Timetree     | [Coalescent multiplication ordering diverges from v0 without empirical validation](M-timetree-coalescent-multiplication-ordering-diverges-from-v0.md) |
| Medium     | Timetree     | [Date column header matching breaks on hash](M-timetree-date-header-hash.md)                                                                        |
| Medium     | Timetree     | [Internal node dates missing in nexus for input branch length mode](M-timetree-internal-dates-missing-input-bl.md)                                  |
| Medium     | Timetree     | [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md)                                                                                |
| Medium     | Timetree     | [Positional likelihood metric differs from v0](M-timetree-positional-likelihood-metric.md)                                                          |
| Medium     | Timetree     | [Skyline coalescent uses Nelder-Mead instead of SLSQP](M-timetree-skyline-nelder-mead-optimizer.md)                                                 |
| Medium     | Timetree     | [Skyline coalescent timing may diverge from v0](M-timetree-skyline-timing-v0-divergence.md)                                                         |
| Medium     | Timetree     | [Marginal dense golden master node key mismatch on ebola_20](M-timetree-dense-golden-master-node-mismatch.md)                                       |
| Medium     | Timetree     | [Marginal dense timetree inference disproportionately slow for mpox dataset](M-timetree-marginal-dense-mpox-slow.md)                                |
| Medium     | Timetree     | [Marginal timetree node times can violate topology](M-timetree-marginal-node-times-can-violate-topology.md)                                         |
| Medium     | Timetree     | [Sparse timetree convergence tracking compares empty sequences](M-timetree-sparse-convergence-empty-sequences.md)                                   |
| Negligible | Ancestral    | [Ancestral command does not produce auspice JSON](N-ancestral-auspice-json-not-produced.md)                                                         |
| Negligible | Ancestral    | [Sparse marginal passes still use remove/insert pattern](N-ancestral-sparse-remove-insert-pattern.md)                                               |
| Negligible | Core         | [Zero branch length clamping](N-core-branch-length-clamping.md)                                                                                     |
| Negligible | Dates        | [Imprecise date upper bound not capped at present](N-dates-imprecise-upper-bound-not-capped.md)                                                     |
| Negligible | Distribution | [Formula discretization errors silently swallowed](N-distribution-formula-silent-discretization.md)                                                 |
| Negligible | I/O          | [Time-based branch lengths not written to Newick/Nexus output](N-io-time-based-branch-lengths-not-implemented.md)                                   |
| Negligible | I/O          | [shared graph writer missing PhyloXML, Auspice, and UShER MAT formats](N-io-write-graph-files-missing-formats.md)                                   |
| Negligible | I/O          | [Large datasets require all sequences in memory simultaneously](N-io-large-dataset-memory-constraint.md)                                            |
| Negligible | I/O          | [Multi-segment genome input not wired](N-io-multi-segment-genome-input.md)                                                                          |
| Negligible | I/O          | [Newick writer defaults to 3 significant digits, truncating branch lengths](N-io-nwk-writer-3-sigfig-default-truncates-precision.md)                |
| Negligible | I/O          | [Graphviz DOT writer uses Newick precision defaults for edge labels](N-io-graphviz-dot-branch-length-precision.md)                                  |
| Negligible | I/O          | [Name reconciliation logic duplicated across subsystems](N-io-name-reconciliation-duplicated.md)                                                    |
| Negligible | I/O          | [Edge annotations from Newick not wired into EdgeFromNwk](N-io-edge-annotation-wiring-not-implemented.md)                                           |
| Negligible | Optimize     | [update_marginal traverses graph twice for mixed partitions](N-optimize-double-graph-traversal-update-marginal.md)                                  |
| Negligible | Optimize     | [initial_guess_mixed allocates Vec\<Sub\> per edge for count only](N-optimize-initial-guess-alloc.md)                                               |
| Negligible | Optimize     | [Dense/sparse equivalence test bounds undocumented](N-optimize-equivalence-bounds-undocumented.md)                                                  |
| Negligible | Optimize     | [Optimize command accepts only a single alignment](N-optimize-multi-alignment-input.md)                                                             |
| Negligible | Optimize     | [Dense optimize iteration is slow](N-optimize-dense-iteration-slow.md)                                                                              |
| Negligible | Optimize     | [Multi-modal surface counterexample with $\ell'(0) < 0$ is not constructed](N-optimize-multimodal-counterexample-unreproduced.md)                   |
| Negligible | Optimize     | [Grid search resolution (100 points) is unverified](N-optimize-grid-search-resolution-unverified.md)                                                |
| Negligible | Optimize     | [Grid search upper bound caps at 0.5 subs/site](N-optimize-grid-search-upper-bound-capped.md)                                                       |
| Negligible | Optimize     | [`reconcile_zero_boundary` grid extent argument is unverified](N-optimize-reconcile-grid-extent-unverified.md)                                      |
| Negligible | Optimize     | [Optimize pre-step uses default Brent tolerance instead of v0's coarse tolerance](N-optimize-pre-step-coarse-tolerance.md)                          |
| Negligible | Optimize     | [Topology cleanup in optimize loop uses Fitch substitutions instead of marginal MAP substitutions](N-optimize-topology-cleanup-fitch-vs-ml-subs.md) |
| Negligible | Timetree     | [--dates not required, misleading error when omitted](N-timetree-dates-not-required.md)                                                             |
| Negligible | Timetree     | [Dead CLI flags in timetree](N-timetree-dead-cli-flags.md)                                                                                          |
| Negligible | Timetree     | [gtr.json missing for --branch-length-mode=input](N-timetree-gtr-json-missing-input-bl.md)                                                          |
| Negligible | Timetree     | [timetree.json missing coalescent and skyline parameters](N-timetree-json-missing-coalescent.md)                                                    |
| Negligible | Timetree     | [Missing output files compared to v0](N-timetree-missing-output-files.md)                                                                           |
| Negligible | Timetree     | [Missing skyline output files](N-timetree-missing-skyline-output.md)                                                                                |
| Negligible | Timetree     | [--n-branches-posterior returns error](N-timetree-n-branches-posterior-unimplemented.md)                                                            |
| Negligible | Timetree     | [Negative coalescent Tc accepted without validation](N-timetree-negative-coalescent-tc.md)                                                          |
| Negligible | Timetree     | [write_node_dates() is a todo!() stub](N-timetree-node-dates-output-unimplemented.md)                                                               |
| Negligible | Timetree     | [--plot-rtt and --plot-tree return error](N-timetree-plot-unimplemented.md)                                                                         |
| Negligible | Timetree     | [--keep-polytomies and --resolve-polytomies no conflicts_with declaration](N-timetree-polytomy-flags-no-conflict.md)                                |
| Negligible | Timetree     | [Stochastic polytomy resolution not implemented](N-timetree-stochastic-polytomy-unimplemented.md)                                                   |
| Negligible | Timetree     | [Auspice JSON output missing mutations, branch confidence, and genome annotations](N-timetree-auspice-json-incomplete.md)                           |
| Negligible | Timetree     | [Node data date string can differ by one day at year-fraction boundaries](N-timetree-node-data-date-string-fp-boundary.md)                          |
| Negligible | Timetree     | [Root node omits placeholder branch-length fields in node data JSON](N-timetree-node-data-root-branch-fields-omitted.md)                            |
| Negligible | Timetree     | [Rerooted root and polytomy-resolution nodes stay unnamed](N-timetree-unnamed-root-after-reroot.md)                                                 |
| Negligible | Timetree     | [Polytomy resolution numerical robustness](N-timetree-polytomy-numerical-robustness.md)                                                             |
| Negligible | Timetree     | [Polytomy resolution test improvements](N-timetree-polytomy-test-improvements.md)                                                                   |
| Medium     | Timetree     | [Branch distribution grid uses uniform spacing](M-timetree-branch-grid-uniform-resolution.md)                                                       |
| Medium     | Timetree     | [Timetree inference in input mode collapses internal-node dates to Empty](M-timetree-inference-input-mode-date-collapse.md)                         |
| Medium     | Timetree     | [sum_coalescent_cost silently clamps negative branch lengths](M-timetree-coalescent-branch-length-clamp.md)                                         |
| Medium     | Clock        | [Clock command discards eleven CLI arguments without wiring](M-clock-dead-cli-arguments.md)                                                         |
| Medium     | Clock        | [Clock command collapses date intervals to their mean](M-clock-date-interval-collapsed-to-mean.md)                                                  |
| Medium     | Discrete     | [Discrete partition uniform_profile(0) produces inf-filled profile](M-discrete-missing-zero-states-inf.md)                                          |
| Medium     | Mugration    | [Mugration optimize_gtr_rate restores mu without restoring profiles](M-mugration-gtr-rate-restore-inconsistency.md)                                 |
| Medium     | GTR          | [Dense GTR golden master bakes v1 root-state filtering into the v0 oracle](M-gtr-dense-root-filter-golden-master-self-validating.md)                |
| Negligible | Mugration    | [Mugration count_transitions lacks direct unit test](N-mugration-count-transitions-untested.md)                                                     |
| Negligible | Timetree     | [Branch grid extent uses base clock_rate, not effective rate](N-timetree-branch-grid-gamma-omitted.md)                                              |
| Negligible | Marginal     | [Marginal forward pass zero-divisor floor converts structural zeros](N-marginal-forward-zero-divisor-floor.md)                                      |
| Negligible | Timetree     | [Convergence metric silently excludes failed coalescent likelihood](N-timetree-convergence-metric-excludes-coalescent.md)                           |
| Medium     | Ancestral    | [Sparse variant skips initialize_marginal before update_marginal](M-ancestral-sparse-missing-initialize-marginal.md)                                |
| Medium     | Marginal     | [normalize_inplace NEG_INFINITY masks all contributions](M-marginal-normalize-neg-infinity-masks-total.md)                                          |
| Medium     | Timetree     | [Timetree convergence metric deficiencies](M-timetree-convergence-metric-deficiencies.md)                                                           |
| Medium     | Timetree     | [Timetree confidence interval computation deficiencies](M-timetree-confidence-interval-deficiencies.md)                                             |
| Medium     | Coalescent   | [Coalescent leaf survival formula sign convention requires investigation](M-coalescent-leaf-survival-sign-convention.md)                            |
| Medium     | Coalescent   | [Coalescent subsystem time notation conflict](M-coalescent-time-notation-conflict.md)                                                               |
| Medium     | Inference    | [Inference forward/backward pass asymmetries](M-inference-forward-backward-asymmetry.md)                                                            |
| Medium     | CLI          | [CLI help text defects, inconsistencies, and UX violations](M-cli-help-text-defects.md)                                                             |
| Medium     | Core         | [Remaining architectural debt after domain module extraction](M-core-remaining-architectural-debt-after-extraction.md)                              |
| Negligible | Timetree     | [Timetree indel rate recomputed per refinement iteration](N-timetree-indel-rate-recomputed-per-iteration.md)                                        |
| Negligible | Numerical    | [Numerical stability magic constants](N-numerical-stability-magic-constants.md)                                                                     |
| Negligible | Coalescent   | [Coalescent skyline robustness deficiencies](N-coalescent-skyline-robustness.md)                                                                    |
| Negligible | Coalescent   | [Coalescent contribution evaluation allocates for scalar merger rates](N-coalescent-contribution-evaluation-allocation-overhead.md)                 |
| Negligible | Core         | [Error suppression via unwrap_or_default and silent fallbacks](N-error-suppression-unwrap-or-defaults.md)                                           |
| Negligible | Alphabet     | [Alphabet serialization format design](N-alphabet-serialization-format-design.md)                                                                   |
| Negligible | Core         | [Code quality and convention violations](N-code-quality-conventions.md)                                                                             |
| Negligible | Core         | [Architectural debt (documented, low priority)](N-architectural-debt-documented.md)                                                                 |
| Negligible | Core         | [Production unwrap/expect/assert audit](N-production-unwrap-expect-audit.md)                                                                        |
| Negligible | Test         | [Test coverage gaps across production functions](N-test-coverage-gaps.md)                                                                           |
| Negligible | Test         | [Test quality deficiencies](N-test-quality-deficiencies.md)                                                                                         |
| Negligible | Test         | [Tests with unnecessary filesystem dependency](N-test-filesystem-dependency.md)                                                                     |
| Negligible | Test         | [Loose tolerances in test_gaussian_product.rs](N-test-gaussian-product-loose-tolerances.md)                                                         |
| Negligible | GTR          | [GTR site-specific interpolation tolerance requires investigation](N-gtr-site-specific-interpolation-tolerance.md)                                  |
| Negligible | Timetree     | [Coalescent integration test uses grossly loose tolerance](N-timetree-coalescent-integration-grossly-loose-tolerance.md)                            |
| Negligible | Core         | [EdgeToGraphviz trait has four identical implementations](N-core-edge-to-graphviz-identical-impls.md)                                               |
| Negligible | Validation   | [ValidationRunner print methods duplicated across three runner implementations](N-validation-runner-print-method-boilerplate.md)                    |
| Negligible | Validation   | [TestCase trait field accessors duplicated across five test suites](N-validation-test-case-accessor-boilerplate.md)                                 |
| Negligible | Graph        | [Parallel traversal partition write locks serialize the frontier](N-graph-parallel-traversal-partition-lock-contention.md)                          |

## Cross-references

- [Design documents](../_raw/) - specifications (work items come from here too, not just v0)
- [Unimplemented algorithms](../algo/unimplemented.md) - algorithm details for unimplemented features
- [Feature inventory](../features/README.md) - v0/v1 feature parity tracking
- [Proposals](../proposals/README.md) - new features pre-implementation
