# Coalescent backward pass grid explosion

Timetree with coalescent prior stalls due to grid size explosion in `convolution_function_function` during the backward propagation in iteration 2+. Severity is dataset-dependent:

- `zika/20 --coalescent=10.0`: hangs (no output after 25+ minutes, disabled in smoke tests)
- `ebola/20 --coalescent=10.0`: 2m41s wall for a 20-tip dataset (completes but slow)
- `ebola/20 --coalescent-skyline`: 5m16s wall for a 20-tip dataset (completes but slow)
- `flu/h3n2/200 --coalescent=10.0`: 1m11s wall for a 200-tip dataset (completes but slow)
- `flu/h3n2/20 --coalescent=10.0`: 2s (unaffected)
- `--coalescent-opt` runs are unaffected (1s) because they use a different Tc estimation path

## Repro

```bash
# Hangs (default --max-iter=2):
./dev/docker/run .build/docker/release/treetime timetree \
  --tree=data/zika/20/tree.nwk --dates=data/zika/20/metadata.tsv \
  --outdir=tmp/timetree/zika/20/coalescent data/zika/20/aln.fasta.xz \
  --coalescent=10.0

# Completes (~2s):
./dev/docker/run .build/docker/release/treetime timetree \
  --tree=data/zika/20/tree.nwk --dates=data/zika/20/metadata.tsv \
  --outdir=tmp/timetree/zika/20/coalescent data/zika/20/aln.fasta.xz \
  --coalescent=10.0 --max-iter=1
```

## Mechanism

1. Coalescent contributions are `Distribution::Formula` with wide range (~265 years for zika/20).

2. `multiply_formula_function()` in [`packages/treetime-distribution/src/distribution_ops/multiply.rs`](../../packages/treetime-distribution/src/distribution_ops/multiply.rs) discretizes the Formula on the Function's grid, using `n_points = b.len()` over the overlap range. When a narrow Function (range ~2 years, 1000 points, dx ~0.002) is multiplied with the wide Formula, the result inherits dx ~0.002.

3. `convolution_function_function()` in [`packages/treetime-distribution/src/distribution_ops/convolve.rs`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs) resamples both inputs to `dx = min(dx_a, dx_b)`. The branch distribution (range ~200 years, dx ~0.2) is resampled from 1000 to ~100,000 points to match the fine child dx.

4. Output size is `a_len + b_len - 1`, reaching 100K-567K points per convolution.

5. This compounds across tree depth. Each level inherits fine dx from children's multiplied distributions.

### Diagnostic evidence (iteration 2, release, zika/20, `--coalescent=10.0`)

Convolution output sizes grow per tree level:

```
76K, 82K, 85K, 100K, 103K, 110K, 114K, 567K
```

Deepest node:

```
a: len=1000, range=[2013.44, 2013.79], dx=3.53e-4
b: len=1000, range=[-200, -0.02], dx=0.2
  resampled b: 566,660 points
  output: 567,659 points
```

### Why iteration 1 works but iteration 2 hangs

- Initial backward pass: leaf time distributions are `Point` (trivial convolution).
- After forward pass + iteration 1: internal nodes get Function distributions with moderate grids.
- Iteration 2: coalescent Formula \* iteration-1 Function (already somewhat fine) produces very fine dx, which cascades through the tree.

## Fix direction

Cap the intermediate grid size in `convolution_function_function`. The result is already resampled back to `max(dx_a, dx_b)` (line 180-181 in convolve.rs), so capping the intermediate dx preserves accuracy while preventing explosion. A max resample size of ~10K points per operand would keep convolutions fast.

## Related issues

- [Branch distribution grid uses uniform spacing](M-timetree-branch-grid-uniform-resolution.md) - same underlying grid architecture
- [Marginal dense timetree inference disproportionately slow for mpox dataset](M-timetree-marginal-dense-mpox-slow.md) - same convolution resampling mechanism with slow-evolving organisms
