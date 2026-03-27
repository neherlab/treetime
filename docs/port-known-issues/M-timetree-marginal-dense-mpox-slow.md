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

`compute_branch_length_distribution()` in [`packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs) creates a 200-point branch-length grid via `create_simple_grid()`, then converts to the time domain by dividing by `effective_clock_rate = clock_rate * gamma`. Different edges have different gamma values, producing time-domain grids with different spacings.

`convolution_function_function()` in [`packages/treetime-distribution/src/distribution_ops/convolve.rs`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs) resamples both inputs to `dx = min(dx_a, dx_b)` before convolving. When one distribution has very fine spacing (long branch, high gamma) and another has coarse spacing (short branch, low gamma), the resampled arrays become large. The convolution output has `len_a + len_b - 1` points, then is resampled back to the coarser grid.

For mpox with its low clock rate, the time-domain grid spacings can differ by large factors across edges, causing repeated expensive resamplings during the backward pass tree traversal.

### v0 handling

v0's `BranchLenInterpolator` ([`packages/legacy/treetime/treetime/branch_len_interpolator.py`](../../packages/legacy/treetime/treetime/branch_len_interpolator.py)) uses adaptive non-uniform grids:

- Short branches (mutation_length < min(1e-5, 0.1 \* one_mutation)): concatenated linear + quadratic grid
- Long branches: logarithmic near zero + quadratic + quadratic tail

v0's `Distribution` class uses FFT-based convolution with adaptive grid merging, which avoids the resampling blow-up by working in a common grid convention.

v1 uses uniform grids throughout (`create_simple_grid` produces `Array1::linspace`), and the convolution resamples to the finer grid. This design is simpler but does not adapt to the branch length regime.

## Paths to investigate

- Profile the backward pass for mpox_clade_ii_20 to confirm the bottleneck is in `convolution_function_function` resampling
- Measure actual grid spacing ratios across edges for mpox vs flu datasets
- Evaluate whether adaptive grid construction (matching v0's Poisson/Gaussian regime switch) would reduce the spacing ratio
- Evaluate whether capping the resampling ratio (e.g. refusing to resample by more than 10x) would preserve accuracy while bounding cost
- Evaluate FFT-based convolution for the Function-Function case (v0 uses FFT via `numpy.fft`)

## Related issues

- [Golden master runner tests missing internal node times for 5 datasets](M-timetree-gm-runner-missing-internal-times.md) - the mpox test also fails with value mismatch, but the 44-minute runtime is a separate concern
- Internal node dates missing at scale (resolved) - may have shared root cause of degenerate distributions for slow-evolving lineages
