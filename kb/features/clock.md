# Clock Inference

## Clock Regression

- [x] Backward regression pass
  - [x] Leaf contribution (date + divergence)
  - [x] Outlier leaf contribution path
  - [x] Internal node aggregation
  - [x] Root assignment
- [x] Forward regression pass for reroot search
- [x] Ordinary least squares with unit tip weights (v0 splits a tip's unit variance along a terminal branch; see [kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))
- [ ] Covariation mode (`--covariation`, `--tip-slack`, variance parameters): removed in v1 ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))
- [ ] Rate and root-date standard deviations in `clock` output (v0: `--covariation`): see [kb/issues/M-clock-command-reports-no-rate-uncertainty.md](../issues/M-clock-command-reports-no-rate-uncertainty.md)

## Rerooting

- [x] Exact least-squares root position on every edge (v0: grid and bounded Brent per edge; see [kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))
- [x] Reroot rejects candidates with non-positive inferred clock rate
- [x] `--keep-root` (disable rerooting; v0 still reroots in its clock filter, see [kb/v0-errata/clock-keep-root-ignored-by-clock-filter.md](../v0-errata/clock-keep-root-ignored-by-clock-filter.md))

## Reroot Modes (4 modes)

- [x] Least-squares (minimize RTT regression residuals)
- [x] Min-dev (fixed-zero-rate least squares; minimize RTT variance independently of sampling-date values; see [kb/v0-errata/clock-min-dev-fixed-slope-score.md](../v0-errata/clock-min-dev-fixed-slope-score.md))
- [x] Oldest (root at oldest node)
- [x] Tip/MRCA mode (`--outgroup`)

## Edge Split Optimization

- [x] Closed-form optimum of the quadratic score on each edge
- [ ] Configurable 1D optimizers (`--branch-split-*`): removed in v1, the optimum is exact ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))

## Edge Operations for Rerooting

- [x] Snap to source endpoint (split fraction below $10^{-5}$, as v0)
- [x] Snap to target endpoint (split fraction above $1 - 10^{-5}$, as v0)
- [x] Edge split (create new node at split point)
- [x] Invert every edge on path to new root
- [x] Edge merge (remove trivial nodes, remove old root when it becomes trivial)

## Clock Filtering

- [x] Outlier prefilter when `--clock-filter > 0`
- [x] Divergence assignment (root-to-tip accumulation)
- [x] IQD calculation (Q3-Q1)
- [x] Outlier marking (threshold \* IQD)
- [ ] Local clock filter method (`--clock-filter-method=local` in v0)
- [ ] `--prune-outliers` (v0 clock command removes outliers from output tree)

## Clock Model Output

- [x] Rate, intercept, chisq, R-value
- [ ] Regression covariance: not reported; the regression treats tips as independent ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))
- [x] Rerooted Newick tree (plain node names; v0 fuses input support into them, see [kb/v0-errata/clock-newick-support-fused-into-node-names.md](../v0-errata/clock-newick-support-fused-into-node-names.md))
- [x] Graph JSON and Graphviz DOT
- [x] Clock model JSON
- [x] Regression CSV (populated for all nodes, includes predicted date, clock deviation, outlier state)
- [x] SVG, PNG, and PDF RTT charts
- [x] Terminal chart printing on TTY

## CLI Options

- [x] `--clock-rate` (fixed rate bypasses estimation)
- [x] `--allow-negative-rate` (for midpoint rooting)
- [ ] `--tip-slack`, `--covariation`, `--sequence-length`: removed in v1 ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))
- [x] `--reroot` / `--outgroup` / `--keep-root`
- [ ] `--aln`, `--gtr`, `--gtr-params`, `--vcf-reference`, `--branch-length-mode`, `--method-anc`, `--prune-short`: removed in v1, because the regression uses the branch lengths of the input tree and v0 reads the alignment only for `--covariation` ([kb/decisions/clock-ordinary-least-squares-regression.md](../decisions/clock-ordinary-least-squares-regression.md))
- [ ] `--seed` (not accepted: clock has no random step; see [kb/decisions/cli-seed-on-commands-with-random-steps.md](../decisions/cli-seed-on-commands-with-random-steps.md))
- [ ] Tree inference from alignment (help text mentions it, not implemented)
