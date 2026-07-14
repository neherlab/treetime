# Marginal dense timetree inference disproportionately slow for mpox dataset

The `test_gm_runner_marginal_dense` golden master test for `mpox_clade_ii_20` (20 tips) takes 2647 seconds (44 minutes) before failing with a value mismatch. For comparison, `flu_h3n2_20` (also 20 tips) completes in under 10 seconds. The 250x slowdown is disproportionate to the tree size and requires investigation.

## Repro

```bash
# Fast (flu, 20 tips): ~7 seconds
./dev/docker/run cargo test -p treetime -- test_gm_runner_marginal_dense::case_2_flu_h3n2_20 --test-threads=1

# Slow (mpox, 20 tips): ~44 minutes
./dev/docker/run cargo test -p treetime -- test_gm_runner_marginal_dense::case_4_mpox_clade_ii_20 --test-threads=1
```

The full CLI also shows the disparity:

```bash
# Fast
./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/flu/h3n2/20/tree.nwk --dates=data/flu/h3n2/20/metadata.tsv \
  --outdir=tmp/timetree/flu/h3n2/20 data/flu/h3n2/20/aln.fasta.xz

# Slow
./dev/docker/run ./dev/dev r treetime -- timetree \
  --tree=data/mpox/clade-ii/20/tree.nwk --dates=data/mpox/clade-ii/20/metadata.tsv \
  --outdir=tmp/timetree/mpox/clade-ii/20 data/mpox/clade-ii/20/aln.fasta.xz
```

## Scientific context

Mpox (monkeypox virus, MPXV) is a double-stranded DNA virus with a ~197 kb genome and a substitution rate around 1-2 x 10^-5 substitutions/site/year. This is roughly 200x slower than influenza H3N2 (RNA virus, ~3-5 x 10^-3 subs/site/year). The low substitution rate means:

- Branch lengths in the tree are very small (few substitutions per site)
- Many branches have zero or near-zero observed mutations
- The branch length probability distributions are sharply peaked near zero (Poisson regime)
- Clock rate estimates are small, producing wide time-domain grids when converting `time = branch_length / (clock_rate * gamma)`

Similar slow-mutation-rate organisms that would trigger the same issue: tuberculosis (Mycobacterium tuberculosis, ~1 x 10^-7 subs/site/year), other DNA viruses, and bacterial pathogens with long generation times.

## Suspected mechanism

The performance bottleneck is in the backward pass (`propagate_distributions_backward`) where child time distributions are convolved with negated branch distributions, then multiplied across children.

### Grid spacing ratio blow-up

`compute_branch_length_distribution()` in [`packages/treetime/src/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/timetree/inference/branch_length_likelihood.rs) creates a 200-point branch-length grid via `create_simple_grid()`, then converts to the time domain by dividing by `effective_clock_rate = clock_rate * gamma`. Different edges have different gamma values, producing time-domain grids with different spacings.

`convolution_function_function()` in [`packages/treetime-distribution/src/distribution_ops/convolve.rs`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs) resamples both inputs to `dx = min(dx_a, dx_b)` before convolving. When one distribution has very fine spacing (long branch, high gamma) and another has coarse spacing (short branch, low gamma), the resampled arrays become large. The convolution output has `len_a + len_b - 1` points, then is resampled back to the coarser grid.

For mpox with its low clock rate, the time-domain grid spacings can differ by large factors across edges, causing repeated expensive resamplings during the backward pass tree traversal.

### Compounding factors

Two properties of mpox compound the grid blow-up:

- Long genome (197 kb). Dense marginal reconstruction processes every alignment column at every node. The per-iteration cost scales as O(nodes _ L _ n_states^2). With L = 197,000 vs 1,800 for flu, each marginal pass is ~110x more expensive before any grid effects.
- Long genome + low rate = tiny `one_mutation`. `one_mutation = 1/L` is ~5e-6 for mpox vs ~5.5e-4 for flu. The grid minimum (`one_mutation * 0.1`) and grid density scale with this value, producing very fine grid spacing in subs/site that then maps to a very wide time-domain grid when divided by the low clock rate.

The net effect: mpox edges produce time-domain distributions with both fine spacing (from the small `one_mutation`) and wide extent (from `MAX_BRANCH_TIME / clock_rate`). When convolved with a child distribution from a different edge with different gamma, the spacing ratio can exceed 100x, inflating resampled arrays to tens of thousands of points.

### v1 convolution: O(n^2) with uniform grids

v1's `convolution_function_function` (`convolve.rs:136`) follows this sequence:

1. Resample both inputs to `dx = min(dx_a, dx_b)` via linear interpolation
2. Direct convolution via `ndarray-conv` (`ConvMode::Full`) - O(n_a \* n_b) where n_a, n_b are resampled lengths
3. Resample output back to `max(dx_a, dx_b)`

When `dx_a / dx_b` is large (grid spacing ratio blow-up), the finer-grid input grows proportionally, making the O(n^2) convolution expensive. The temporary output array (`n_a + n_b - 1` points at the fine spacing) can also be large.

### v0 handling: FFT convolution with adaptive grids and delta shortcuts

v0's `NodeInterpolator.convolve_fft` (`node_interpolator.py:161`) uses three strategies that avoid the blow-up:

1. Delta shortcut. When the node distribution is much narrower than the branch distribution (ratio < 1/FFT_FWHM_GRID_SIZE), treat the node as a delta function and shift the branch distribution. This avoids convolution entirely for sharply peaked distributions (common for well-constrained leaf dates).

2. FFT convolution. For the general case, both distributions are evaluated on a common uniform grid with spacing `dt = max(one_mutation * 0.005, min(fwhm_node, fwhm_branch) / fft_grid_size)`. The convolution is computed as `irfft(rfft(f) * rfft(g))` at O(n log n) cost. The grid is sized at `2 * max(support_a, support_b) / dt`, adapting to the actual support width.

3. Adaptive non-uniform branch grids. `BranchLenInterpolator` (`branch_len_interpolator.py:37`) builds grids adapted to the branch length regime:
   - Short branches (mutation_length < min(1e-5, 0.1 \* one_mutation)): linear + quadratic grid concentrated near zero
   - Long branches: logarithmic near zero (5 points spanning 20 orders of magnitude), quadratic left of peak, quadratic right to 3\*sigma, quadratic far tail to MAX_BRANCH_LENGTH
   - The grid is non-uniform, placing resolution where the likelihood varies and using sparse coverage in the tails

v1 uses uniform `Array1::linspace` grids throughout (`create_simple_grid` produces 200 uniform points from `one_mutation * 0.1` to `max(3 * center, MAX_BRANCH_TIME * clock_rate)`). This design is simpler but does not adapt to the branch length regime and forces the convolution to resample between grids with very different spacings.

## Paths to investigate

- Profile the backward pass for mpox_clade_ii_20 to confirm the bottleneck is in `convolution_function_function` resampling
- Measure actual grid spacing ratios across edges for mpox vs flu datasets
- Evaluate whether adaptive grid construction (matching v0's Poisson/Gaussian regime switch) would reduce the spacing ratio
- Evaluate whether capping the resampling ratio (e.g. refusing to resample by more than 10x) would preserve accuracy while bounding cost
- Evaluate FFT-based convolution for the Function-Function case (v0 uses FFT via `numpy.fft`)

## Related issues

- Source: [kb/issues/M-timetree-marginal-dense-mpox-slow.md](../issues/M-timetree-marginal-dense-mpox-slow.md) -- delete after full resolution
- [Golden master runner tests missing internal node times for 5 datasets](../issues/M-timetree-gm-runner-missing-internal-times.md) -- the mpox test also fails with value mismatch, but the 44-minute runtime is a separate concern
- The same slow-evolving regime can also amplify degenerate or weakly informative time distributions, so this performance issue is a plausible co-factor when internal-node inference fails on low-rate datasets
