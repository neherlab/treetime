# `--no-indels` flag to disable indel contributions to branch-length optimization

## Feature

Add a `--no-indels` CLI flag for `optimize` and `timetree` that disables the Poisson indel contribution to the branch-length likelihood. Default: indels enabled (current behavior). When set: the optimizer uses substitution-only likelihood, matching v0 and standard phylogenetic tools (RAxML, IQ-TREE, PhyML, BEAST).

The indel contribution is a v1 design-doc feature ([optimize-indel-contribution-to-likelihood](../port-intentional-changes/optimize-indel-contribution-to-likelihood.md)). The toggle enables:

- v0 parity testing: golden master comparison requires matching v0's substitution-only likelihood
- Controlled experiments: isolating the effect of indel contributions on branch length estimates
- Datasets with alignment artifacts where indel counts are unreliable

## Code paths to gate

### Optimize loop (`run_optimize_loop`)

| Location                                              | Role                                                                             |
| ----------------------------------------------------- | -------------------------------------------------------------------------------- |
| `packages/treetime/src/commands/optimize/run.rs#L397` | `estimate_indel_rate()` called per iteration in `compute_iteration_likelihood()` |
| `packages/treetime/src/commands/optimize/run.rs#L398` | `total_indel_log_lh()` summed into total likelihood                              |
| `packages/treetime/src/commands/optimize/run.rs#L359` | `run_optimize_mixed_with_indel_rate()` passes rate to per-edge optimizer         |

### Per-edge optimization (`run_optimize_mixed_with_indel_rate`)

| Location                                                           | Role                                                                  |
| ------------------------------------------------------------------ | --------------------------------------------------------------------- |
| `packages/treetime/src/commands/optimize/optimize_unified.rs#L552` | Per-edge loop queries `edge_indel_count()` and passes to optimizer    |
| `packages/treetime/src/commands/optimize/optimize_unified.rs#L247` | `evaluate_with_indels()` adds Poisson term to substitution likelihood |
| `packages/treetime/src/commands/optimize/method_newton.rs`         | Newton methods receive `indel_count` and `indel_rate`                 |
| `packages/treetime/src/commands/optimize/method_brent.rs`          | Brent methods receive `indel_count` and `indel_rate`                  |

### Timetree branch distributions

| Location                                                                        | Role                                                                           |
| ------------------------------------------------------------------------------- | ------------------------------------------------------------------------------ |
| `packages/treetime/src/commands/timetree/inference/runner.rs#L93`               | `estimate_indel_rate()` for branch distribution grid                           |
| `packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs` | `compute_branch_length_distribution()` receives `indel_count` and `indel_rate` |

### CLI args

| Location                                          | Role                                       |
| ------------------------------------------------- | ------------------------------------------ |
| `packages/treetime/src/commands/optimize/args.rs` | `TreetimeOptimizeArgs` has no indel toggle |

## Implementation sketch

Add a `--no-indels` boolean flag (default false) to `TreetimeOptimizeArgs`. When set:

1. In `run_optimize_loop`: skip `estimate_indel_rate()` (use `0.0`), skip `total_indel_log_lh()` (use `0.0`)
2. In `run_optimize_mixed_with_indel_rate`: skip `edge_indel_count()` queries (use `0`), pass `indel_rate = 0.0`
3. In timetree `compute_branch_distributions_marginal_mode`: skip `estimate_indel_rate()` (use `0.0`), skip `edge_indel_count()` (use `0`)

When both `indel_count = 0` and `indel_rate = 0.0`, `poisson_indel_log_lh()` returns `OptimizationMetrics::default()` (all zeros), cleanly dropping the Poisson term from the likelihood, its derivatives, and the Newton/Brent optimizers.

The flag should also be added to the timetree command args (timetree reuses the same branch-length optimization infrastructure).
