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

| Severity   | Scope        | Issue                                                                                                                                                              |
| ---------- | ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| High       | App          | [Application transport contracts diverge across clients](H-app-transport-contracts-diverge-across-clients.md)                                                      |
| High       | App          | [Application UI ignores inputs and displays synthetic scientific results](H-app-ui-displays-synthetic-results-and-ignores-inputs.md)                               |
| High       | Benchmark    | [Benchmark compared-revision execution trust is undecided](H-benchmark-compared-revision-execution-trust-undecided.md)                                             |
| High       | Benchmark    | [Benchmark paths cross shell, filesystem, and Markdown boundaries unsafely](H-benchmark-paths-cross-shell-filesystem-and-markdown-boundaries.md)                   |
| High       | Clock        | [ClockSet dateless-leaf contribution biases regression](H-clock-clockset-dateless-leaf-bias.md)                                                                    |
| High       | Core         | [Command modules contain shared operations that belong in domain layers](H-core-command-module-shared-ops-entanglement.md)                                         |
| High       | Core         | [Library crate is not consumable by non-CLI clients](H-core-multi-client-architecture-library-purity.md)                                                           |
| High       | Distribution | [Result-returning distribution operations panic on Formula operands](H-distribution-result-api-panics-on-formula.md)                                               |
| High       | Graph        | [Indexed payload extraction exposes invalid shared graph state](H-graph-indexed-payload-extraction-exposes-invalid-state.md)                                       |
| High       | Homoplasy    | [Homoplasy command is unimplemented](H-homoplasy-command-unimplemented.md)                                                                                         |
| High       | I/O          | [Auspice v2 output omits required updated metadata](H-io-auspice-v2-required-updated-missing.md)                                                                   |
| High       | I/O          | [UShER MAT export fabricates global reference nucleotides from branch-local alleles](H-io-usher-ref-nuc-uses-parent-allele.md)                                     |
| High       | Marginal     | [Marginal cavity sentinel loses impossible-factor multiplicity](H-marginal-cavity-sentinel-loses-impossible-factor-multiplicity.md)                                |
| High       | Timetree     | [Coalescent model may be built before node times are recomputed after a topology change](H-timetree-coalescent-events-incomplete-after-topology-change.md)         |
| High       | Timetree     | [No tree inference from alignment](H-timetree-tree-inference-unimplemented.md)                                                                                     |
| Medium     | Ancestral    | [Dense-sparse log-likelihood divergence](M-ancestral-dense-sparse-divergence.md)                                                                                   |
| Medium     | Ancestral    | [Fitch recurrence is not minimum parsimony on multifurcations](M-ancestral-fitch-polytomy-recurrence-not-minimum.md)                                               |
| Medium     | Ancestral    | [GenBank annotation format not supported for AA reconstruction](M-ancestral-genbank-annotation-unsupported.md)                                                     |
| Medium     | Ancestral    | [Marginal reconstruction uses plain probability space](M-ancestral-marginal-probability-space.md)                                                                  |
| Medium     | Ancestral    | [Missing initialize_marginal before update_marginal in sparse and AA dense paths](M-ancestral-sparse-missing-initialize-marginal.md)                               |
| Medium     | Ancestral    | [Sparse root invariance violation](M-ancestral-sparse-root-invariance.md)                                                                                          |
| Medium     | Ancestral    | [Sparse variable-site alphabet mismatch](M-ancestral-sparse-alphabet-mismatch.md)                                                                                  |
| Medium     | App          | [Browser file ownership is undecided](M-app-browser-file-ownership-undecided.md)                                                                                   |
| Medium     | App          | [Frontend verification excludes the renderer](M-app-frontend-verification-excludes-renderer.md)                                                                    |
| Medium     | Benchmark    | [Benchmark reports mix revisions and are not reproducible](M-benchmark-reports-mix-revisions-and-are-not-reproducible.md)                                          |
| Medium     | CLI          | [CLI help text defects, inconsistencies, and UX violations](M-cli-help-text-defects.md)                                                                            |
| Medium     | Clock        | [Clock command collapses date intervals to their mean](M-clock-date-interval-collapsed-to-mean.md)                                                                 |
| Medium     | Clock        | [Clock command discards eight CLI arguments without wiring](M-clock-dead-cli-arguments.md)                                                                         |
| Medium     | Clock        | [Clock covariation overdispersion hardcoded](M-clock-covariation-overdispersion.md)                                                                                |
| Medium     | Clock        | [Clock filter residual computation differs from v0](M-clock-filter-residual-parity.md)                                                                             |
| Medium     | Clock        | [MinDev reroot uses wrong objective (EstimatedRate instead of FixedRate(0))](M-clock-mindev-wrong-objective.md)                                                    |
| Medium     | Clock        | [Parallel clock CSV row order is nondeterministic](M-clock-parallel-output-order-nondeterministic.md)                                                              |
| Medium     | Coalescent   | [Coalescent bad_branch handling is inconsistent across code paths](M-coalescent-bad-branch-handling-inconsistent.md)                                               |
| Medium     | Coalescent   | [Coalescent edge collection bypasses NaN dates and retains unreachable multiplicity fallback](M-coalescent-edge-collection-nan-bypass-and-unreachable-fallback.md) |
| Medium     | Coalescent   | [Coalescent weight applied to Range distributions evaluates only at endpoints](M-coalescent-range-endpoint-only-interpolation.md)                                  |
| Medium     | Core         | [Mutation representation and format projection are inconsistent](M-core-mutation-representation-and-format-projection-inconsistent.md)                             |
| Medium     | Core         | [Partition creation is hardcoded per-command instead of configured](M-core-partition-init-orchestration-duplication.md)                                            |
| Medium     | Core         | [Remaining architectural debt after domain module extraction](M-core-remaining-architectural-debt-after-extraction.md)                                             |
| Medium     | Dates        | [Column auto-detection gaps in CSV readers](M-dates-column-auto-detection-gaps.md)                                                                                 |
| Medium     | Discrete     | [Discrete partition constructor accepts zero states](M-discrete-missing-zero-states-inf.md)                                                                        |
| Medium     | Distribution | [Distribution normalization erases formula and grid errors](M-distribution-normalization-erases-errors.md)                                                         |
| Medium     | Distribution | [Distribution product grid resolution diverges from v0](M-distribution-product-grid-resolution-diverges-from-v0.md)                                                |
| Medium     | Distribution | [Distribution Range accepts invalid supports and amplitudes](M-distribution-range-unchecked-domain.md)                                                             |
| Medium     | Distribution | [Distribution support-boundary semantics are unresolved](M-distribution-support-boundary-semantics-unresolved.md)                                                  |
| Medium     | Distribution | [Mixed-support distribution operations lose exact intersection boundaries](M-distribution-mixed-support-operations-lose-exact-intersection-boundaries.md)          |
| Medium     | Distribution | [Neg-log likely time selects the least-likely ordinate](M-distribution-neglog-likely-time-selects-maximum.md)                                                      |
| Medium     | Distribution | [Plain distribution division applies an unscaled fixed divisor floor](M-distribution-plain-division-fixed-floor.md)                                                |
| Medium     | Graph        | [Reroot edge splitting lacks validation and failure atomicity](M-graph-reroot-split-validation-and-atomicity.md)                                                   |
| Medium     | Graph        | [Target topology order silently collapses duplicate labels](M-topology-order-duplicate-label-collapse.md)                                                          |
| Medium     | Grid         | [Grid constructors accept non-finite and non-representable spacing](M-grid-invalid-numeric-domain.md)                                                              |
| Medium     | GTR          | [Dense GTR golden master bakes v1 root-state filtering into the "v0" oracle](M-gtr-dense-root-filter-golden-master-self-validating.md)                             |
| Medium     | GTR          | [Per-site rate variation is not inferred or wired end to end](M-gtr-per-site-rate-variation.md)                                                                    |
| Medium     | I/O          | [Auspice entropy projection perturbs the Shannon definition](M-io-auspice-entropy-perturbs-shannon-definition.md)                                                  |
| Medium     | I/O          | [Auspice trait fields are silently dropped when malformed](M-io-auspice-trait-confidence-silently-coerced.md)                                                      |
| Medium     | I/O          | [GFF3 CDS phase is ignored and not validated](M-io-gff3-cds-phase-ignored.md)                                                                                      |
| Medium     | I/O          | [PhyloXML boolean properties are silently coerced](M-io-phyloxml-booleans-silently-coerced.md)                                                                     |
| Medium     | I/O          | [Sequence attachment has $O(n^2)$ complexity](M-io-sequence-attachment-quadratic.md)                                                                               |
| Medium     | I/O          | [Sequence-to-node name matching is unreliable](M-io-sequence-name-matching-unreliable.md)                                                                          |
| Medium     | I/O          | [Tree-backed output order is inconsistent](M-io-tree-backed-output-order-inconsistent.md)                                                                          |
| Medium     | I/O          | [UShER input collapses explicit zero branch length into absence](M-io-usher-zero-branch-length-collapsed.md)                                                       |
| Medium     | I/O          | [UShER MAT output drops unrepresentable mutations implicitly](M-io-usher-mat-mutation-loss-is-implicit.md)                                                         |
| Medium     | I/O          | [VCF input and output are unimplemented](M-io-vcf-input-output-unimplemented.md)                                                                                   |
| Medium     | Inference    | [Fallible parallel inference passes partially commit state](M-inference-fallible-parallel-passes-partially-commit.md)                                              |
| Medium     | Inference    | [Inference forward/backward pass asymmetries](M-inference-forward-backward-asymmetry.md)                                                                           |
| Medium     | Marginal     | [normalize_inplace NEG_INFINITY masks all contributions](M-marginal-normalize-neg-infinity-masks-total.md)                                                         |
| Medium     | Mugration    | [Mugration rejects metadata entries not present in tree](M-mugration-extra-metadata-names-rejected.md)                                                             |
| Medium     | Mugration    | [Mugration golden master parity with v0](M-mugration-iterative-gtr.md)                                                                                             |
| Medium     | Mugration    | [Mugration optimize_gtr_rate restores mu without restoring profiles](M-mugration-gtr-rate-restore-inconsistency.md)                                                |
| Medium     | Optimize     | [Bifurcating root optimization: independent per-edge vs joint arc optimization](M-optimize-root-bifurcating-independent-vs-joint.md)                               |
| Medium     | Optimize     | [Indel event count reduced after composition on merged edges](M-optimize-indel-event-count-after-composition.md)                                                   |
| Medium     | Optimize     | [Optimize omits v0 post-convergence short-branch pruning](M-optimize-short-branch-pruning-unimplemented.md)                                                        |
| Medium     | Optimize     | [Optimize per-branch lengths diverge from v0 fixture beyond 10% relative tolerance](M-optimize-gm-per-branch-divergence.md)                                        |
| Medium     | Optimize     | [Optimizer and marginal propagation use incompatible branch-length scales when the GTR rate scalar differs from one](M-optimize-gtr-mu-coordinate-mismatch.md)     |
| Medium     | Sparse       | [Sparse marginal fixed_counts zero-multiplicity produces NaN](M-sparse-marginal-zero-multiplicity-nan.md)                                                          |
| Medium     | Timetree     | [--keep-root rejects negative clock rate](M-timetree-keep-root-rejects-negative-clock-rate.md)                                                                     |
| Medium     | Timetree     | [--method-anc ignored in timetree](M-timetree-method-anc-ignored.md)                                                                                               |
| Medium     | Timetree     | [Branch distribution grid uses uniform spacing](M-timetree-branch-grid-uniform-resolution.md)                                                                      |
| Medium     | Timetree     | [Constant-$T_c$ optimization can report success without a parameter](M-timetree-constant-tc-success-without-parameter.md)                                          |
| Medium     | Timetree     | [Date column header matching breaks on hash](M-timetree-date-header-hash.md)                                                                                       |
| Medium     | Timetree     | [Golden master runner tests missing internal node times for 5 datasets](M-timetree-gm-runner-missing-internal-times.md)                                            |
| Medium     | Timetree     | [Internal node dates missing in nexus for input branch length mode](M-timetree-internal-dates-missing-input-bl.md)                                                 |
| Medium     | Timetree     | [Marginal dense golden master node key mismatch on ebola_20](M-timetree-dense-golden-master-node-mismatch.md)                                                      |
| Medium     | Timetree     | [Marginal dense timetree inference disproportionately slow for mpox dataset](M-timetree-marginal-dense-mpox-slow.md)                                               |
| Medium     | Timetree     | [Marginal timetree node times can violate topology](M-timetree-marginal-node-times-can-violate-topology.md)                                                        |
| Medium     | Timetree     | [Positional likelihood metric differs from v0](M-timetree-positional-likelihood-metric.md)                                                                         |
| Medium     | Timetree     | [Skyline optimization aborts on node-time inversions instead of degrading gracefully](M-timetree-skyline-aborts-on-node-time-inversion.md)                         |
| Medium     | Timetree     | [Skyline objective and optimizer diverge from v0](M-timetree-skyline-nelder-mead-optimizer.md)                                                                     |
| Medium     | Timetree     | [Sparse timetree convergence tracking compares empty sequences](M-timetree-sparse-convergence-empty-sequences.md)                                                  |
| Medium     | Timetree     | [Timetree backward pass combines child time messages in plain probability space](M-timetree-backward-pass-plain-space-underflow.md)                                |
| Medium     | Timetree     | [Timetree forward pass skips uncertain leaf posterior refinement](M-timetree-forward-pass-skips-uncertain-leaf-refinement.md)                                      |
| Medium     | Timetree     | [Timetree confidence interval computation deficiencies](M-timetree-confidence-interval-deficiencies.md)                                                            |
| Medium     | Timetree     | [Timetree convergence metric deficiencies](M-timetree-convergence-metric-deficiencies.md)                                                                          |
| Medium     | Timetree     | [Timetree inference in input mode collapses internal-node dates to Empty](M-timetree-inference-input-mode-date-collapse.md)                                        |
| Medium     | Timetree     | [Timetree skyline coalescent timing may diverge from v0](M-timetree-skyline-timing-v0-divergence.md)                                                               |
| Medium     | Timetree     | [TreeTime tree-output inference metadata is incomplete](M-timetree-tree-output-inference-metadata-incomplete.md)                                                   |
| Negligible | Alphabet     | [Alphabet serialization format design](N-alphabet-serialization-format-design.md)                                                                                  |
| Negligible | Ancestral    | [Ancestral Auspice output is incomplete and method-dependent](N-ancestral-auspice-json-not-produced.md)                                                            |
| Negligible | Ancestral    | [Fitch site classification regresses parallel command performance](N-ancestral-fitch-site-classification-parallel-regression.md)                                   |
| Negligible | Ancestral    | [Fitch transmission filtering can produce an empty state set](N-ancestral-fitch-empty-transmission-state-set.md)                                                   |
| Negligible | Ancestral    | [Marginal ndarray kernels allocate avoidable intermediates](N-ancestral-marginal-array-kernels-allocate.md)                                                        |
| Negligible | Ancestral    | [Parallel sparse leaf setup has no defined error atomicity contract](N-ancestral-parallel-sparse-leaf-error-atomicity-unverified.md)                               |
| Negligible | Ancestral    | [Parallel sparse leaf setup may regress single-thread ancestral runtime](N-ancestral-parallel-sparse-leaf-single-thread-regression.md)                             |
| Negligible | Ancestral    | [Parallel sparse leaf setup validation covers only a narrow success path](N-ancestral-parallel-sparse-leaf-validation-coverage.md)                                 |
| Negligible | Ancestral    | [Sparse reconstruction can lose node state and suppress invariant failures](N-ancestral-sparse-remove-insert-pattern.md)                                           |
| Negligible | App          | [Application toolchain rules conflict with package manifests](N-app-toolchain-rules-conflict-with-manifests.md)                                                    |
| Negligible | Array        | [Owned ndarray signatures force projection and propagation copies](N-array-owned-signatures-force-projection-copies.md)                                            |
| Negligible | Build        | [Reusable library selects the native linear-algebra backend](N-build-library-selects-linear-algebra-backend.md)                                                    |
| Negligible | Clock        | [Clock SVG output is nondeterministic for an untraced reason](N-clock-svg-output-nondeterministic-untraced.md)                                                     |
| Negligible | Clock        | [ClockSet::chisq() numerically unstable for near-zero determinant](N-clock-chisq-near-zero-determinant.md)                                                         |
| Negligible | Clock        | [v1 clock regression produces all-negative rates where v0 finds positive](N-clock-regression-all-negative-rate.md)                                                 |
| Negligible | Coalescent   | [Coalescent time-scale coordinate is not type-enforced](N-coalescent-time-scale-coordinate-not-type-enforced.md)                                                   |
| Negligible | Coalescent   | [Skyline constant extrapolation policy is unapproved](N-coalescent-skyline-extrapolation-policy-undecided.md)                                                      |
| Negligible | Coalescent   | [Skyline grid construction accepts invalid endpoint arrays](N-coalescent-skyline-grid-validation-incomplete.md)                                                    |
| Negligible | Coalescent   | [Skyline merger-rate quadrature lacks an accuracy contract](N-coalescent-skyline-quadrature-contract-undecided.md)                                                 |
| Negligible | Coalescent   | [Skyline simplex initialization is scale-independent and one-sided](N-coalescent-skyline-simplex-initialization-undecided.md)                                      |
| Negligible | Core         | [Architectural debt (documented, low priority)](N-architectural-debt-documented.md)                                                                                |
| Negligible | Core         | [Code quality and convention violations](N-code-quality-conventions.md)                                                                                            |
| Negligible | Core         | [Error suppression via unwrap_or_default and silent fallbacks](N-error-suppression-unwrap-or-defaults.md)                                                          |
| Negligible | Core         | [Production unwrap/expect/assert audit](N-production-unwrap-expect-audit.md)                                                                                       |
| Negligible | Core         | [Zero branch length clamping](N-core-branch-length-clamping.md)                                                                                                    |
| Negligible | Dates        | [Imprecise date upper bound not capped at present](N-dates-imprecise-upper-bound-not-capped.md)                                                                    |
| Negligible | Distribution | [Formula discretization errors silently swallowed](N-distribution-formula-silent-discretization.md)                                                                |
| Negligible | Distribution | [Mixed-NaN distribution semantics are undefined](N-distribution-mixed-nan-policy-undecided.md)                                                                     |
| Negligible | Docs         | [Auto-partitioning report has citation metadata and notation defects](N-kb-auto-partitioning-citation-metadata-and-notation.md)                                    |
| Negligible | Docs         | [Documentation references and source claims are not verifiable](N-doc-reference-and-source-integrity.md)                                                           |
| Negligible | Docs         | [Knowledge base contains internal links with missing targets](N-kb-internal-links-have-missing-targets.md)                                                         |
| Negligible | Docs         | [Mathematical notation is inconsistent and leaves symbols undeclared](N-doc-mathematical-notation-inconsistent.md)                                                 |
| Negligible | Docs         | [Unimplemented algorithm inventory has incomplete citation structure](N-kb-unimplemented-algorithm-citation-structure.md)                                          |
| Negligible | Graph        | [common_ancestor missing direct tests](N-common-ancestor-missing-direct-tests.md)                                                                                  |
| Negligible | Graph        | [Dependency-frontier schedulers are duplicated](N-graph-dependency-frontier-schedulers-duplicated.md)                                                              |
| Negligible | Graph        | [iter_children_arc scans full node list per visited node](N-graph-traverse-quadratic-iter-children-arc.md)                                                         |
| Negligible | GTR          | [GTR site-specific interpolation tolerance requires investigation](N-gtr-site-specific-interpolation-tolerance.md)                                                 |
| Negligible | GTR          | [Site-specific GTR inference lacks end-to-end test from real tree data](N-gtr-site-specific-e2e-inference-test.md)                                                 |
| Negligible | GTR          | [Site-specific GTR not integrated into partition system](N-gtr-site-specific-partition-integration.md)                                                             |
| Negligible | I/O          | [Auspice number formatting lacks an independent Augur golden master](N-io-auspice-number-formatting-missing-augur-golden-master.md)                                |
| Negligible | I/O          | [Edge annotations from Newick not wired into EdgeFromNwk](N-io-edge-annotation-wiring-not-implemented.md)                                                          |
| Negligible | I/O          | [Graphviz DOT writer uses Newick precision defaults for edge labels](N-io-graphviz-dot-branch-length-precision.md)                                                 |
| Negligible | I/O          | [Large datasets require all sequences in memory simultaneously](N-io-large-dataset-memory-constraint.md)                                                           |
| Negligible | I/O          | [Multi-segment genome input not wired](N-io-multi-segment-genome-input.md)                                                                                         |
| Negligible | I/O          | [Name reconciliation logic duplicated across subsystems](N-io-name-reconciliation-duplicated.md)                                                                   |
| Negligible | I/O          | [Newick writer defaults to 3 significant digits, truncating branch lengths](N-io-nwk-writer-3-sigfig-default-truncates-precision.md)                               |
| Negligible | I/O          | [Nexus parser handles a limited subset](N-io-nexus-parser-limited-subset.md)                                                                                       |
| Negligible | I/O          | [PhyloXML mutation property contract is undecided](N-io-phyloxml-mutation-property-contract-undecided.md)                                                          |
| Negligible | I/O          | [Root branch length silently discarded during Newick parsing](N-io-newick-root-branch-length-discarded.md)                                                         |
| Negligible | I/O          | [Shared TreeIR architecture has not received an explicit design decision](N-io-tree-ir-architecture-unapproved.md)                                                 |
| Negligible | I/O          | [Time-based branch lengths not written to Newick/Nexus output](N-io-time-based-branch-lengths-not-implemented.md)                                                  |
| Negligible | Marginal     | [Marginal forward pass zero-divisor floor converts structural zeros](N-marginal-forward-zero-divisor-floor.md)                                                     |
| Negligible | Marginal     | [Sparse cavity-profile compression lacks an error contract](N-sparse-marginal-cavity-compression-unverified.md)                                                    |
| Negligible | Mugration    | [Mugration transition counting lacks a complete discrete hand oracle](N-mugration-count-transitions-untested.md)                                                   |
| Negligible | Numerical    | [Numerical stability magic constants](N-numerical-stability-magic-constants.md)                                                                                    |
| Negligible | Optimize     | [`reconcile_zero_boundary` grid extent argument is unverified](N-optimize-reconcile-grid-extent-unverified.md)                                                     |
| Negligible | Optimize     | [Dense optimize iteration is slow](N-optimize-dense-iteration-slow.md)                                                                                             |
| Negligible | Optimize     | [Dense/sparse equivalence test bounds undocumented](N-optimize-equivalence-bounds-undocumented.md)                                                                 |
| Negligible | Optimize     | [Grid search resolution (100 points) is unverified](N-optimize-grid-search-resolution-unverified.md)                                                               |
| Negligible | Optimize     | [Grid search upper bound caps at 0.5 subs/site](N-optimize-grid-search-upper-bound-capped.md)                                                                      |
| Negligible | Optimize     | [initial_guess_mixed allocates Vec<Sub> per edge for count only](N-optimize-initial-guess-alloc.md)                                                                |
| Negligible | Optimize     | [Multi-modal surface counterexample with a negative likelihood derivative at zero is not constructed](N-optimize-multimodal-counterexample-unreproduced.md)        |
| Negligible | Optimize     | [Optimize command accepts only a single alignment](N-optimize-multi-alignment-input.md)                                                                            |
| Negligible | Optimize     | [Optimize pre-step uses default Brent tolerance instead of v0's coarse tolerance](N-optimize-pre-step-coarse-tolerance.md)                                         |
| Negligible | Optimize     | [Sparse branch optimization state parity lacks a coefficient oracle](N-optimize-sparse-state-parity-unverified.md)                                                 |
| Negligible | Optimize     | [Topology cleanup in optimize loop uses Fitch substitutions instead of marginal MAP substitutions](N-optimize-topology-cleanup-fitch-vs-ml-subs.md)                |
| Negligible | Optimize     | [update_marginal traverses the graph twice for mixed partitions](N-optimize-double-graph-traversal-update-marginal.md)                                             |
| Negligible | Partition    | [`infer_dense()` stub always returns false](N-representation-infer-dense-stub.md)                                                                                  |
| Negligible | Partition    | [Dense and sparse partition types have structural and naming asymmetries](N-representation-dense-sparse-partition-asymmetry.md)                                    |
| Negligible | Reroot       | [Clock rerooting duplicates the generic root-search implementation](N-reroot-clock-search-duplicates-generic-module.md)                                            |
| Negligible | Reroot       | [Leaf variance offset convention diverges from v0 on terminal edges](N-reroot-leaf-variance-offset-diverges-from-v0.md)                                            |
| Negligible | Reroot       | [No v0 golden master test for min-dev reroot](N-reroot-missing-v0-golden-master.md)                                                                                |
| Negligible | Reroot       | [Reroot edge-stats maps use BTreeMap where Vec indexing fits](N-reroot-btreemap-edge-stats-perf.md)                                                                |
| Negligible | Reroot       | [Reroot split-optimizer default diverges from v0](N-reroot-split-optimizer-default-diverges-from-v0.md)                                                            |
| Negligible | Reroot       | [resolve_tip_keys error and ambiguity paths untested](N-reroot-tip-resolution-untested-errors.md)                                                                  |
| Negligible | Reroot       | [Tip-name resolution duplicated between optimize and clock reroot](N-reroot-duplicated-tip-name-resolution.md)                                                     |
| Negligible | Test         | [Floating-point assertion tolerance is hidden by a variable](N-test-floating-tolerance-hidden-by-variable.md)                                                      |
| Negligible | Test         | [Loose tolerances in test_gaussian_product.rs](N-test-gaussian-product-loose-tolerances.md)                                                                        |
| Negligible | Test         | [Test coverage gaps across production functions](N-test-coverage-gaps.md)                                                                                          |
| Negligible | Test         | [Test quality deficiencies](N-test-quality-deficiencies.md)                                                                                                        |
| Negligible | Test         | [Unit tests depend on repository-level production datasets](N-test-filesystem-dependency.md)                                                                       |
| Negligible | Timetree     | [--dates not required, misleading error when omitted](N-timetree-dates-not-required.md)                                                                            |
| Negligible | Timetree     | [--keep-polytomies and --resolve-polytomies no conflicts_with declaration](N-timetree-polytomy-flags-no-conflict.md)                                               |
| Negligible | Timetree     | [--n-branches-posterior returns error](N-timetree-n-branches-posterior-unimplemented.md)                                                                           |
| Negligible | Timetree     | [--plot-rtt and --plot-tree return error](N-timetree-plot-unimplemented.md)                                                                                        |
| Negligible | Timetree     | [Branch grid extent uses base clock_rate, not effective rate](N-timetree-branch-grid-gamma-omitted.md)                                                             |
| Negligible | Timetree     | [Convergence metric silently excludes failed coalescent likelihood](N-timetree-convergence-metric-excludes-coalescent.md)                                          |
| Negligible | Timetree     | [Dead CLI flags in timetree](N-timetree-dead-cli-flags.md)                                                                                                         |
| Negligible | Timetree     | [gtr.json missing for --branch-length-mode=input](N-timetree-gtr-json-missing-input-bl.md)                                                                         |
| Negligible | Timetree     | [Missing output files compared to v0](N-timetree-missing-output-files.md)                                                                                          |
| Negligible | Timetree     | [Missing skyline output files](N-timetree-missing-skyline-output.md)                                                                                               |
| Negligible | Timetree     | [Negative coalescent Tc accepted without validation](N-timetree-negative-coalescent-tc.md)                                                                         |
| Negligible | Timetree     | [Node data `date` string can differ by one day at year-fraction boundaries](N-timetree-node-data-date-string-fp-boundary.md)                                       |
| Negligible | Timetree     | [Polytomy resolution numerical robustness](N-timetree-polytomy-numerical-robustness.md)                                                                            |
| Negligible | Timetree     | [Polytomy resolution test improvements](N-timetree-polytomy-test-improvements.md)                                                                                  |
| Negligible | Timetree     | [Rerooted root and polytomy-resolution nodes stay unnamed](N-timetree-unnamed-root-after-reroot.md)                                                                |
| Negligible | Timetree     | [Root node omits placeholder branch-length fields in node data JSON](N-timetree-node-data-root-branch-fields-omitted.md)                                           |
| Negligible | Timetree     | [Stochastic polytomy event generator not implemented](N-timetree-stochastic-polytomy-unimplemented.md)                                                             |
| Negligible | Timetree     | [Timetree indel rate recomputed per refinement iteration](N-timetree-indel-rate-recomputed-per-iteration.md)                                                       |
| Negligible | Timetree     | [timetree.json missing coalescent and skyline parameters](N-timetree-json-missing-coalescent.md)                                                                   |
| Negligible | Timetree     | [Timetree --vary-rate is unimplemented](N-timetree-vary-rate-unimplemented.md)                                                                                     |
| Negligible | Timetree     | [write_node_dates() is a todo!() stub](N-timetree-node-dates-output-unimplemented.md)                                                                              |

## Cross-references

- [Design documents](../_raw/) - specifications (work items come from here too, not just v0)
- [Unimplemented algorithms](../algo/unimplemented.md) - algorithm details for unimplemented features
- [Feature inventory](../features/README.md) - v0/v1 feature parity tracking
- [Proposals](../proposals/README.md) - new features pre-implementation
