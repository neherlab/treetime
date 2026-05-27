# Implement root profile sampling for ancestral reconstruction

Add stochastic sampling from the posterior profile at the root node during sequence assignment. Augur always uses this (`sample_from_profile='root'`). v1 uses deterministic argmax only.

## Changes

1. Define `SampleMode` enum: `Argmax`, `SampleRoot`, `SampleAll`.

2. Modify `reconstruct_map_seq` in `packages/treetime/src/partition/marginal_sparse.rs` to accept `SampleMode` and an optional `&mut impl Rng`. At root when `SampleRoot`: cumsum the profile row, generate uniform random, find first index where cumsum >= random.

3. Same change for the dense path in `packages/treetime/src/partition/marginal_dense.rs`.

4. Thread the RNG from `--seed` CLI arg through `ancestral_reconstruction_marginal` in `packages/treetime/src/ancestral/marginal.rs`.

5. Wire `SampleMode::SampleRoot` as default for the ancestral command (matching augur's behavior). Add CLI flag `--sample-from-profile` with values `none`, `root`, `all`.

6. Default to `SampleRoot` to match augur. `--seed` controls reproducibility.

## Related issues

- Source: [M-ancestral-root-sampling-not-implemented.md](../issues/M-ancestral-root-sampling-not-implemented.md)
