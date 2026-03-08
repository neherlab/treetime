# Zero-length branches cause panic

Zero-length or near-zero branches in the input tree propagate through belief
propagation, producing distribution grids with non-uniform spacing. The
interpolation layer panics with `x array must be uniformly spaced`. Affected
datasets: dengue_20, lassa_L_20, mpox_clade_ii_20, rsv_a_20, tb_20. These
datasets are commented out in the timetree runner tests:

- [test_gm_runner_poisson.rs#L19](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_poisson.rs#L19)
- [test_gm_runner_marginal_dense.rs#L31](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_dense.rs#L31)
- [test_gm_runner_marginal_sparse.rs#L32](../../packages/treetime/src/commands/timetree/inference/__tests__/test_gm_runner/test_gm_runner_marginal_sparse.rs#L32)

## Locations

- Assertion: [`grid.rs#L83`](../../packages/treetime-grid/src/grid.rs#L83)
- Convolution output: [`convolve.rs#L133`](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L133)
- Call site: [`backward_pass.rs#L65`](../../packages/treetime/src/commands/timetree/inference/backward_pass.rs#L65)

## Workaround

Affected datasets can produce partial timetree output using
`--branch-length-mode=input`, bypassing the grid construction code. The output
is tips-only (see
[Internal node dates missing for coalescent modes](M-timetree-internal-dates-missing-coalescent.md))
but the clock rate may be valid (e.g. 0.000380 for TB).

## Repro

```bash
# tb
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 \
  --tree=data/tb/20/tree.nwk \
  --dates=data/tb/20/metadata.tsv \
  --name-column=strain --date-column=date \
  --outdir=tmp/repro-grid-tb data/tb/20/aln.fasta.xz
# x array must be uniformly spaced

# lassa
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 \
  --tree=data/lassa/L/20/tree.nwk \
  --dates=data/lassa/L/20/metadata.tsv \
  --outdir=tmp/repro-grid-lassa data/lassa/L/20/aln.fasta.xz
# Clock rate 5.399e-4 valid but low mutation rate causes grid artifacts

# mpox
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 \
  --tree=data/mpox/clade-ii/20/tree.nwk \
  --dates=data/mpox/clade-ii/20/metadata.tsv \
  --outdir=tmp/repro-grid-mpox data/mpox/clade-ii/20/aln.fasta.xz
# Clock rate 6.4e-5 (very low)
```

## Related issues

- [Zero branch length clamping](N-core-branch-length-clamping.md) related but
  occurs at a different point in the pipeline - the clamping in
  `fix_branch_length` applies during ancestral reconstruction, not during
  timetree distribution construction
- [Negative clock rate causes assertion failure](H-timetree-negative-clock-rate-assertion.md)
  blocks the same datasets (rsv, dengue) through a different mechanism
