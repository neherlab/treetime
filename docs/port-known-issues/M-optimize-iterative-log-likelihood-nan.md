# Optimize iterative log-likelihood becomes `-inf`/NaN on larger datasets

The optimize command completes the initial sparse-message pass without crashing, but later iterations still drive the reported log-likelihood to `-inf`/NaN on larger datasets. The failure appears in both dense and sparse modes, so the remaining defect is no longer the old sparse-only `combine_messages()` underflow.

## Scope

- Reproduced on `ebola/100`
- Also observed on `flu/h3n2/200`
- Affects both `--dense=false` and dense mode

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- optimize --dense=false \
  --tree=data/ebola/100/tree.nwk --aln=data/ebola/100/aln.fasta.xz \
  --outdir=tmp/optimize/sparse/ebola/100
```

## Root cause

The forward-pass branch update in [packages/treetime/src/representation/partition/marginal_passes.rs](../../packages/treetime/src/representation/partition/marginal_passes.rs) divides the parent profile by the child message:

- [packages/treetime/src/representation/partition/marginal_passes.rs](../../packages/treetime/src/representation/partition/marginal_passes.rs#L235-L268)

When the divisor contains zeros or a numerically degenerate distribution, `numerator / divisor` produces invalid values. The subsequent normalization step takes `norm.ln()` and reuses `dis / norm`, so zero or non-finite `norm` propagates `-inf`/NaN through later iterations.

## Proposed fix

- Guard the forward-pass division against zero and non-finite divisors before `numerator / divisor`
- Define the fallback semantics for degenerate rows, matching the stabilized backward-pass normalization rules
- Add a regression test that runs multiple optimize iterations on a dataset large enough to hit the failure

## Related issues

- [Dense normalize_inplace produces NaN for all-zero probability rows](N-dense-normalize-inplace-zero-row.md)
