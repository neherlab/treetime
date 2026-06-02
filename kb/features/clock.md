# Clock Inference

## Clock Regression

- [x] Backward regression pass
  - [x] Leaf contribution (date + divergence)
  - [x] Outlier leaf contribution path
  - [x] Internal node aggregation
  - [x] Root assignment
- [x] Forward regression pass for reroot search
- [x] Variance propagation (branch length proportional + offset)
- [x] Covariation mode (`--covariation` for regression)
  - [x] Uses `--sequence-length`
  - [x] Uses `--tip-slack` or default `3.0`
  - [x] Builds custom variance parameters
  - [/] Missing `--sequence-length` does not raise explicit error

## Rerooting

- [x] Node scan for best positive-rate root candidate
- [x] Parent-branch split optimization around best node
- [x] Child-branch split optimization around best node
- [x] Reroot rejects candidates with non-positive inferred clock rate
- [x] `--keep-root` (disable rerooting)

## Reroot Modes (4 modes)

- [x] Least-squares (minimize RTT regression residuals)
- [x] Min-dev (minimize RTT variance)
- [x] Oldest (root at oldest node)
- [x] Tip/MRCA mode (`--reroot-tips`)

## Edge Split Optimization (3 methods)

- [x] Grid search (configurable grid points)
- [x] Brent's method (configurable tolerance/max-iter)
- [x] Golden section search (configurable tolerance/max-iter)

## Edge Operations for Rerooting

- [x] Snap to source endpoint
- [x] Snap to target endpoint
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

- [x] Rate, intercept, chisq, R-value, hessian, covariance
- [x] Rerooted Newick tree
- [x] Graph JSON and Graphviz DOT
- [x] Clock model JSON
- [x] Regression CSV (populated for all nodes, includes predicted date, clock deviation, outlier state)
- [x] SVG and PNG RTT charts
- [x] Terminal chart printing on TTY

## CLI Options

- [x] `--clock-rate` (fixed rate bypasses estimation)
- [x] `--allow-negative-rate` (for midpoint rooting)
- [x] `--tip-slack` (terminal node excess variance)
- [x] `--reroot` / `--reroot-tips`
- [ ] `--aln` (parsed but not wired)
- [ ] `--vcf-reference` (parsed but not wired)
- [ ] `--gtr` (parsed but not wired)
- [ ] `--gtr-params` (parsed but not wired)
- [ ] `--branch-length-mode` (parsed but not wired)
- [ ] `--method-anc` (parsed but not wired)
- [ ] `--prune-short` (parsed but not wired)
- [ ] `--seed` (parsed but not wired)
- [ ] Tree inference from alignment (help text mentions it, not implemented)
