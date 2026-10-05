# Knowledge base contains internal links with missing targets

A check of every relative Markdown link target in `kb/` (excluding `kb/_raw/`) finds links whose target file does not exist. The missing targets prevent readers and automated checks from following the KB's evidence graph. Most of them point to source files that moved during refactors.

## Removed `packages/treetime/src/commands/` tree

Command code moved to `packages/app-commands/src/commands/` and the output code to `packages/app-output/src/`. Links to the old `packages/treetime/src/commands/...` paths are in:

- `kb/algo/ancestral.md`
- `kb/algo/optimization.md`
- `kb/algo/unimplemented.md`
- `kb/decisions/ancestral-iterative-gtr-refinement.md`
- `kb/decisions/coalescent-output-schema.md`
- `kb/decisions/command-optimize-standalone.md`
- `kb/decisions/command-prune-standalone.md`
- `kb/decisions/datetime-numeric-date-naming.md`
- `kb/decisions/prune-merge-jukes-cantor-branch-length.md`
- `kb/decisions/timetree-no-zero-branch-collapse-in-loop.md`
- `kb/issues/H-homoplasy-command-unimplemented.md`
- `kb/issues/M-command-output-ownership-is-scattered.md`
- `kb/issues/M-core-mutation-representation-and-format-projection-inconsistent.md`
- `kb/issues/M-output-module-mixes-topology-ordering.md`
- `kb/issues/M-timetree-clock-filter-default-differs-from-v0.md`
- `kb/issues/N-amino-acid-mutation-indel-representation-undecided.md`
- `kb/issues/N-ancestral-auspice-json-not-produced.md`
- `kb/issues/N-datetime-date-and-range-representation-inconsistent.md`
- `kb/issues/N-optimize-topology-cleanup-fitch-vs-ml-subs.md`
- `kb/issues/N-timetree-augur-root-branch-field-tests-stale.md`
- `kb/proposals/amino-acid-mutation-output-substitution-only.md`
- `kb/proposals/optimize-convergence-and-robustness.md`
- `kb/reports/iterative-tree-refinement/8-initial-estimation.md`
- `kb/reports/optimization-methods/4-outer-loops.md`
- `kb/reports/optimization-methods/7-audit.md`
- `kb/reports/optimize-dense-sparse-architecture.md`
- `kb/reports/sparse-subs-accessors.md`
- `kb/reports/zero-branch-length-optimization.md`
- `kb/v0-errata/optimize-signed-convergence-check.md`

## Removed flat `packages/treetime/src/partition/` modules

The partition code moved into the `fitch/`, `marginal/`, `optimize/`, and `storage/` subdirectories. Links to the old flat files (`marginal_core.rs`, `marginal_dense.rs`, `marginal_discrete.rs`, `marginal_helpers.rs`, `marginal_passes.rs`, `marginal_sparse.rs`, `optimize_dense.rs`, `optimize_sparse.rs`, `optimization_contribution.rs`, `dense.rs`, `sparse.rs`, `fitch.rs`, `timetree.rs`, `traits.rs`, `indexed_pass.rs`) are in:

- `kb/algo/ancestral.md`
- `kb/algo/indel-models.md`
- `kb/algo/optimization.md`
- `kb/decisions/ancestral-iterative-gtr-refinement.md`
- `kb/decisions/command-optimize-standalone.md`
- `kb/decisions/mugration-root-state-filtering.md`
- `kb/decisions/optimize-dense-initial-guess-hard-count.md`
- `kb/issues/H-marginal-cavity-sentinel-loses-impossible-factor-multiplicity.md`
- `kb/issues/M-discrete-missing-zero-states-inf.md`
- `kb/issues/M-mugration-tree-leaves-missing-from-metadata-rejected.md`
- `kb/issues/M-partition-probability-profile-field-abbreviated.md`
- `kb/issues/M-partition-sequence-state-field-stutters.md`
- `kb/issues/M-timetree-consumers-read-unconstrained-branch-lengths.md`
- `kb/issues/N-ancestral-deleted-position-invariant-enforced-by-scattered-guards.md`
- `kb/issues/N-marginal-forward-zero-divisor-floor.md`
- `kb/issues/N-mugration-shared-core-alphabet-agnosticism-unguarded.md`
- `kb/issues/N-optimize-sparse-state-parity-unverified.md`
- `kb/issues/N-optimize-topology-cleanup-fitch-vs-ml-subs.md`
- `kb/issues/N-representation-dense-sparse-partition-asymmetry.md`
- `kb/issues/N-sparse-marginal-cavity-compression-unverified.md`
- `kb/proposals/input-name-matching-validation.md`
- `kb/proposals/optimize-convergence-and-robustness.md`
- `kb/reports/2026-07-13_perf-report-parallel-scaling/profiling-parallel-sparse-leaf-setup.md`
- `kb/reports/auto-partitioning.md`
- `kb/reports/codon-substitution-models.md`
- `kb/reports/iterative-tree-refinement/3-tree-likelihood.md`
- `kb/reports/iterative-tree-refinement/4-ancestral-reconstruction.md`
- `kb/reports/optimize-dense-sparse-architecture.md`
- `kb/reports/sparse-subs-accessors.md`
- `kb/reports/stochastic.simulation.md`

## Other removed source files

- `kb/decisions/multi-format-tree-io.md`: `packages/treetime-cli/src/convert/` (`args.rs`, `auspice.rs`, `convert.rs`), `packages/phyloxml/src/types.rs`, `packages/treetime-io/src/auspice.rs`, `packages/treetime-io/src/phyloxml.rs`
- `kb/decisions/coalescent-analytic-tc-optimization.md`, `kb/decisions/optimize-newton-raphson-per-edge.md`, `kb/reports/optimization-methods/7-audit.md`: `packages/treetime/src/coalescent/optimize_tc.rs`
- `kb/decisions/mugration-pseudo-count-initial-pi.md`, `kb/issues/M-discrete-missing-zero-states-inf.md`, `kb/issues/M-mugration-analysis-interface-exposes-policy-wiring.md`: `packages/treetime/src/mugration/mugration.rs`
- `kb/decisions/prune-merge-jukes-cantor-branch-length.md`: `packages/treetime/src/gtr/jc_distance/__tests__/test_jc_distance.rs`
- `kb/issues/H-graph-capability-contracts-silently-discard-state.md`: `packages/treetime/src/clock/clock_graph.rs`
- `kb/issues/M-payloads-own-format-adapters.md`, `kb/reports/sparse-subs-accessors.md`: `packages/treetime/src/payload/` (`ancestral.rs`, `timetree.rs`)
- `kb/issues/N-ancestral-auspice-json-not-produced.md`: `packages/treetime-io/src/graph.rs`
- `kb/issues/N-coalescent-api-unused-state-and-inconsistent-types.md`: `packages/treetime/src/coalescent/precomputed.rs`
- `kb/issues/N-io-auspice-trait-reader-contract-unverified.md`: `packages/treetime-io/src/auspice.rs`
- `kb/reports/auto-partitioning.md`: `packages/treetime/src/gtr/gtr_site_specific.rs`, `packages/treetime/src/gtr/infer_gtr/site_specific.rs`
- `kb/reports/dense-openblas-profiling.md`: `dev/dev`. The script no longer exists, so the report's `./dev/dev bp`, `br`, and `t` commands are stale as well
- `kb/issues/N-doc-treetime-readme-module-layout-stale.md`: `packages/treetime/README.md`

## Repository-root paths used as relative links

`kb/proposals/distribution-log-space-and-hard-soft-boundaries-2.md` links to `packages/...` paths without the `../../` prefix, so the links resolve inside `kb/proposals/`. The affected targets are `distribution_core/distribution.rs`, `distribution_ops/convolve.rs`, `distribution_ops/mass_domain.rs`, `distribution_ops/multiply.rs` of `treetime-distribution`, `grid_fn.rs` of `treetime-grid`, and `backward_pass.rs`, `branch_length_likelihood.rs`, `forward_pass.rs` of `packages/treetime/src/timetree/inference/`.

## Missing knowledge-base entries

- `kb/decisions/command-optimize-standalone.md`: `kb/issues/N-timetree-node-data-confidence-not-emitted.md`
- `kb/decisions/timetree-no-zero-branch-collapse-in-loop.md`, `kb/issues/N-optimize-dense-iteration-slow.md`: `kb/reports/command-relationships/README.md`
- `kb/issues/H-homoplasy-command-unimplemented.md`: `kb/features/README.md`
- `kb/issues/M-timetree-clock-filter-default-differs-from-v0.md`: `kb/issues/H-cli-timetree-config-disables-clock-filter.md`
- `kb/issues/N-ancestral-fitch-site-classification-parallel-scaling-unverified.md`, `kb/reports/2026-07-29_perf-report-work-first-dag-traversal/report.md`: `kb/issues/M-benchmark-reports-mix-revisions-and-are-not-reproducible.md`
- `kb/reports/2026-07-29_perf-report-work-first-dag-traversal/report.md`: `kb/issues/H-benchmark-compared-revision-execution-trust-undecided.md`, `kb/issues/M-clock-parallel-output-order-nondeterministic.md`
- `kb/proposals/distribution-log-space-and-hard-soft-boundaries.md`: `kb/decisions/distribution-intersection-grid-resolution.md`, `kb/decisions/timetree-inference-pass-boundary-tails.md`

## Fix

- Trace moved source files with Git history and update links to their current tracked locations and relevant line anchors
- Replace links to deleted issues, decisions, reports, and proposals with their live successor when one exists; otherwise remove the stale statement or rewrite it to stand alone
- Replace links to untracked build artifacts and removed scripts with stable reproduction instructions or tracked evidence
- Preserve the scientific and design meaning of each surrounding paragraph
- Leave `kb/_raw/` unchanged

## Validation

- A repository-wide Markdown target check confirms that every tracked relative link target exists
