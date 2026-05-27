# Ancestral root sampling not implemented

v1 ancestral reconstruction uses `argmax_first()` for sequence assignment at all nodes including root. v0 supports stochastic sampling from the posterior profile at the root via `sample_from_profile='root'`, which augur always passes.

## Mechanism in v0

`seq_utils.py` `prof2seq()`: when `sample_from_prof=True`, samples from the cumulative distribution at each position using numpy RNG. When `'root'`, samples at root only, argmax everywhere else. Inverse CDF sampling: cumsum the profile row, generate uniform random, find first index where cumsum >= random.

This is orthogonal to joint vs marginal reconstruction. Both paths accept the parameter. It is not related to joint reconstruction removal.

## Effect on output

Only affects positions where the root posterior has near-equal probabilities for multiple states. Changes the root sequence, which changes root mutations (and reference sequence if `--root-sequence` not provided). At positions with a clear winner, sampling almost always picks the same state as argmax.

Augur always passes `sample_from_profile='root'` (`augur/ancestral.py:294`, `augur/refine.py:406`). The RNG is seeded from `--seed` CLI arg. Without `--seed`, non-reproducible.

## v1 status

v1 uses `argmax_first()` everywhere (`packages/treetime/src/partition/marginal_sparse.rs`). Deterministic tie-breaking. No sampling infrastructure.

## Implementation scope

Localized to `reconstruct_map_seq` in `marginal_sparse.rs` and the equivalent in `marginal_dense.rs`:

1. Accept a `SampleMode` enum (`Argmax`, `SampleRoot`, `SampleAll`)
2. For `SampleRoot`: at root node only, sample from profile CDF instead of `argmax_first`
3. Thread an RNG (seeded from `--seed`) through reconstruction
4. ~5 lines of sampling: cumsum profile row, generate uniform, find first index where cumsum >= random

## Locations

- `packages/treetime/src/partition/marginal_sparse.rs` - `reconstruct_map_seq` (sparse path)
- `packages/treetime/src/partition/marginal_dense.rs` - dense equivalent
- `packages/treetime/src/ancestral/marginal.rs` - `ancestral_reconstruction_marginal` (caller)
- v0: `packages/legacy/treetime/treetime/seq_utils.py` `prof2seq()`
