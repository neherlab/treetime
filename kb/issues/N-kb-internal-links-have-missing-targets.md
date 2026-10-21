# Knowledge base contains internal links with missing targets

A check of every relative Markdown link target in `kb/` (outside `kb/_raw/`) finds links whose target file does not exist. The missing targets prevent readers and automated checks from following the KB's evidence graph. Most of them point to source files that moved during refactors.

## Removed flat `packages/treetime/src/partition/` modules

The partition code moved into the `fitch/`, `marginal/`, `optimize/`, and `storage/` subdirectories. Links to the old flat files (`marginal_core.rs`, `marginal_dense.rs`, `marginal_discrete.rs`, `marginal_helpers.rs`, `marginal_passes.rs`, `marginal_sparse.rs`, `optimize_dense.rs`, `optimize_sparse.rs`, `optimization_contribution.rs`, `sparse.rs`, `fitch.rs`, `timetree.rs`, `traits.rs`) are in:

- `kb/algo/indel-models.md`
- `kb/decisions/optimize-dense-initial-guess-hard-count.md`
- `kb/proposals/input-name-matching-validation.md`
- `kb/proposals/optimize-convergence-and-robustness.md`
- `kb/reports/auto-partitioning.md`
- `kb/reports/codon-substitution-models.md`
- `kb/reports/iterative-tree-refinement/3-tree-likelihood.md`
- `kb/reports/iterative-tree-refinement/4-ancestral-reconstruction.md`
- `kb/reports/sparse-subs-accessors.md`
- `kb/reports/stochastic.simulation.md`

## Other removed source files

- `kb/decisions/coalescent-analytic-tc-optimization.md`, `kb/reports/optimization-methods/7-audit.md`: `packages/treetime/src/coalescent/optimize_tc.rs`
- `kb/decisions/prune-merge-jukes-cantor-branch-length.md`: `packages/treetime/src/gtr/jc_distance/__tests__/test_jc_distance.rs`
- `kb/reports/sparse-subs-accessors.md`: `packages/treetime/src/payload/ancestral.rs`
- `kb/reports/stochastic.simulation.md`: `packages/treetime/src/coalescent/edge_data.rs`, `packages/treetime/src/coalescent/events.rs`, `packages/treetime/src/coalescent/integration.rs`
- `kb/reports/auto-partitioning.md`: `packages/treetime/src/gtr/gtr_site_specific.rs`, `packages/treetime/src/gtr/infer_gtr/site_specific.rs`
- `kb/reports/dense-openblas-profiling.md`: `dev/dev`. The script no longer exists, so the report's `./dev/dev bp`, `br`, and `t` commands are stale as well

## Repository-root paths used as relative links

`kb/proposals/distribution-log-space-and-hard-soft-boundaries-2.md` links to `packages/...` paths without the `../../` prefix, so the links resolve inside `kb/proposals/`. The affected targets are `distribution_ops/convolve.rs`, `distribution_ops/mass_domain.rs`, `distribution_ops/multiply.rs` of `treetime-distribution`, `grid_fn.rs` of `treetime-grid`, `backward_pass.rs`, `branch_length_likelihood.rs`, `forward_pass.rs` of `packages/treetime/src/timetree/inference/`, and `packages/treetime/src/timetree/max_posterior_region.rs`.

## Missing knowledge-base entries

- `kb/decisions/timetree-no-zero-branch-collapse-in-loop.md`: `kb/reports/command-relationships/README.md`
- `kb/proposals/distribution-log-space-and-hard-soft-boundaries.md`: `kb/decisions/distribution-intersection-grid-resolution.md`, `kb/decisions/timetree-inference-pass-boundary-tails.md`

## Fix

- Trace moved source files with Git history and update links to their current tracked locations and relevant line anchors
- Replace links to deleted issues, decisions, reports, and proposals with their live successor when one exists; otherwise remove the stale statement or rewrite it to stand alone
- Replace links to untracked build artifacts and removed scripts with stable reproduction instructions or tracked evidence
- Preserve the scientific and design meaning of each surrounding paragraph

## Validation

- A repository-wide Markdown target check confirms that every tracked relative link target exists
