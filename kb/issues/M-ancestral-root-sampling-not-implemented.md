# Ancestral root sampling not implemented

v1 ancestral reconstruction uses `argmax_first()` for sequence assignment at all nodes including root. v0 supports stochastic sampling from the posterior profile at the root via `sample_from_profile='root'`, which augur always passes.

## Mechanism in v0

`seq_utils.py` `prof2seq()`: when `sample_from_prof=True`, samples from the cumulative distribution at each position using numpy RNG. When `'root'`, samples at root only, argmax everywhere else. Inverse CDF sampling: cumsum the profile row, generate uniform random, find first index where cumsum >= random.

This is orthogonal to joint vs marginal reconstruction. Both paths accept the parameter. It is not related to joint reconstruction removal.

## Effect on output

Only affects positions where the root posterior has near-equal probabilities for multiple states. Changes the root sequence, which changes root mutations (and reference sequence if `--root-sequence` not provided). At positions with a clear winner, sampling almost always picks the same state as argmax.

Augur always passes `sample_from_profile='root'` (`augur/ancestral.py:294`, `augur/refine.py:406`). The RNG is seeded from `--seed` CLI arg. Without `--seed`, non-reproducible.

## v1 status

Sampling infrastructure implemented:

- `SampleMode` enum (`Argmax`, `Root`, `All`) in `packages/treetime/src/ancestral/sample.rs`
- `sample_from_profile()` using inverse CDF sampling
- `reconstruct_map_seq_sampled()` in `marginal_sparse.rs` accepts optional RNG
- `prof2seq_sampled()` in `marginal_dense.rs` accepts optional RNG

Remaining: wire `SampleMode` through `ancestral_reconstruction_marginal` and the CLI `--sample-from-profile` flag. Default to `SampleMode::Root` to match augur behavior.

## Locations

- `packages/treetime/src/ancestral/sample.rs` - `SampleMode`, `sample_from_profile()`
- `packages/treetime/src/partition/marginal_sparse.rs` - `reconstruct_map_seq_sampled()`
- `packages/treetime/src/partition/marginal_dense.rs` - `prof2seq_sampled()`
- `packages/treetime/src/ancestral/marginal.rs` - `ancestral_reconstruction_marginal` (pending wiring)
- v0: `packages/legacy/treetime/treetime/seq_utils.py` `prof2seq()`
