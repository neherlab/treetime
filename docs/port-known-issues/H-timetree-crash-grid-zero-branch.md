# Timetree crashes on zero-length branches with grid spacing error

The timetree command crashes with "x array must be uniformly spaced" on datasets containing zero-length branches in the input tree. This affects 5 of 10 test datasets with default parameters and 2 additional datasets with `--keep-root`.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/rsv/a/tree.nwk --dates=data/rsv/a/metadata.tsv \
  --outdir=tmp/timetree/rsv/a data/rsv/a/aln.fasta.xz
```

Any of these datasets crash with default timetree: `rsv/a`, `dengue/500`, `lassa/L`, `mpox/clade-ii`, `tb`.

With `--keep-root`: `flu/h3n2/20`, `flu/h3n2/200` also crash.

## Error

```
Error: x array must be uniformly spaced
```

From `Grid::from_array()` at [`packages/treetime-grid/src/grid.rs#L82-L83`](../../packages/treetime-grid/src/grid.rs#L82-L83).

## Root cause

The crash occurs during the timetree backward pass when distributions with different grid spacings are combined via convolution or multiplication.

### Crash path

1. `compute_branch_distributions_marginal_mode()` ([`packages/treetime/src/commands/timetree/inference/runner.rs#L76`](../../packages/treetime/src/commands/timetree/inference/runner.rs#L76)) creates a `Distribution::Function` for each edge via `compute_branch_length_distribution()`.

2. `create_simple_grid()` ([`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L65-L69`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L65-L69)) builds the grid from `one_mutation * 0.1` to `max(center * 3.0, one_mutation * 10.0)`. For zero-length branches (`center = 0.0`), the grid ranges from `one_mutation * 0.1` to `one_mutation * 10.0`. The grid itself is uniform (from `Array1::linspace`).

3. The grid is converted to time domain by dividing by `effective_clock_rate`. Different edges have different `gamma` values, producing grids with different spacings.

4. During `propagate_distributions_backward()` ([`packages/treetime/src/commands/timetree/inference/backward_pass.rs#L17`](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs#L17)), child time distributions are convolved with negated branch distributions (`distribution_convolution`), then combined across children via `distribution_multiplication`.

5. When two `Distribution::Function` values with different grid spacings are convolved or multiplied, downstream `Grid::from_array()` calls receive arrays that are not uniformly spaced and fail the `has_uniform_spacing()` check.

### Why zero-length branches trigger this

Zero-length branches produce branch distributions with extreme shapes (sharp peak near zero branch length). When the branch length grid is converted to the time domain via `time = branch_length / (clock_rate * gamma)`, the resulting time grid spacing differs from that of neighboring branches with non-zero lengths. The combination operations assume compatible grids but receive incompatible ones.

## v0 solution

v0 handles zero-length and short branches through `BranchLenInterpolator` ([`packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102`](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102)):

- Switches between Poisson regime (short branches, `mutation_length <= 0.05`) and Gaussian regime (long branches)
- Uses adaptive grid construction with logarithmic spacing near zero and quadratic tails
- Constructs per-node `NodeInterpolator` distributions that share a common grid convention

v1 uses a single `create_simple_grid()` linear grid for all branches without adaptation for the zero-length case.

## Affected tests

Multiple test files have commented-out cases for the affected datasets:

- [`packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L19-L25`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L19-L25)
- [`packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L32-L38`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L32-L38)
- [`packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L31`](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L31)

## Potential fixes

### F1: Grid compatibility in distribution operations

Make `distribution_multiplication` and `distribution_convolution` resample to a common grid when two `Function` distributions have different spacings. `GridFn::from_arrays_nonuniform()` already handles this case but is not invoked in these code paths.

### F2: Adaptive grid construction (match v0)

Replace `create_simple_grid()` with a `BranchLenInterpolator`-style adaptive grid that uses logarithmic spacing near zero and adapts to branch length regime (Poisson/Gaussian). Ensures all branch distributions share compatible grid conventions.

### F3: Point distribution for zero-length branches

When a branch has zero observed mutations and zero input branch length, use `Distribution::Point(0.0, 1.0)` instead of a grid-based Function distribution. Point distributions bypass grid operations entirely. This would match the input-mode path at [`runner.rs#L168`](../../packages/treetime/src/commands/timetree/inference/runner.rs#L168).

## Related issues

- [Internal node dates missing at scale](M-timetree-internal-dates-missing-scale.md) - shares the root cause of degenerate distributions but affects nodes that survive the grid check
- [Internal node dates missing with bad fixed rate](M-timetree-internal-dates-bad-fixed-rate.md) - wrong clock rate produces similar distribution degeneracy
