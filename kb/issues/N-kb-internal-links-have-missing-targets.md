# Knowledge base contains internal links with missing targets

A check of every relative Markdown link target in `kb/` (outside `kb/_raw/`) finds links whose target file does not exist. The missing targets prevent readers and automated checks from following the KB's evidence graph. Most of them point to source files that moved during refactors.

## Removed flat `packages/treetime/src/partition/` modules

The partition code moved into the `fitch/`, `marginal/`, `optimize/`, and `storage/` subdirectories. Links to the old flat files remain in:

- `kb/reports/iterative-tree-refinement/3-tree-likelihood.md`: `marginal_helpers.rs`, `marginal_passes.rs`
- `kb/reports/iterative-tree-refinement/4-ancestral-reconstruction.md`: `marginal_passes.rs`
- `kb/reports/sparse-subs-accessors.md`: `marginal_passes.rs`, `marginal_sparse.rs`, `optimize_sparse.rs`, `sparse.rs`
- `kb/reports/stochastic.simulation.md`: `marginal_core.rs`, `marginal_discrete.rs`

## Other removed source files

- `kb/reports/optimization-methods/7-audit.md`: `packages/treetime/src/coalescent/optimize_tc.rs`
- `kb/reports/sparse-subs-accessors.md`: `packages/treetime/src/payload/ancestral.rs`
- `kb/reports/stochastic.simulation.md`: `packages/treetime/src/coalescent/edge_data.rs`, `packages/treetime/src/coalescent/events.rs`, `packages/treetime/src/coalescent/integration.rs`
- `kb/issues/N-test-quality-deficiencies.md`: `packages/treetime/src/partition/marginal/sparse/__tests__/test_message.rs`
- `kb/reports/dense-openblas-profiling.md`: `dev/dev`. The script no longer exists, so the report's `./dev/dev bp`, `br`, and `t` commands are stale as well

## Fix

- Trace moved source files with Git history and update links to their current tracked locations and relevant line anchors
- Replace links to removed scripts with stable reproduction instructions or tracked evidence
- Preserve the scientific and design meaning of each surrounding paragraph

## Validation

- A repository-wide Markdown target check confirms that every tracked relative link target exists. The check must parse link destinations in angle brackets as a whole, because DOI URLs such as `<https://doi.org/10.1016/0304-4149(82)...>` contain parentheses
