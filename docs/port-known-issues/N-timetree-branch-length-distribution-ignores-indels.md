# Timetree branch length distribution ignores indels

The timetree command's branch length distribution calculation at [packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L50](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L50) calls `evaluate_mixed_log_lh_only()` directly, which evaluates substitution likelihood only. The Poisson indel contribution added to `run_optimize_mixed()` does not reach this code path.

## Current behavior

`compute_branch_length_distribution()` builds a grid of branch lengths and evaluates the substitution log-likelihood at each point via `evaluate_mixed_log_lh_only()`. This produces a branch length distribution used for time inference. Indel events on the edge do not influence this distribution.

## Expected behavior

The branch length distribution should include the indel contribution, consistent with how `run_optimize_mixed()` evaluates edges. An edge with indels but no substitutions should have a distribution peaked away from zero.

## Impact

Low for typical viral datasets (indels are rare). Visible on branches where the only evolutionary signal is an indel event: the timetree distribution will peak at zero while the optimize command correctly assigns positive length.

## Fix

Pass the indel count and rate to `compute_branch_length_distribution()` and use `evaluate_with_indels_log_lh_only()` (or equivalent) instead of `evaluate_mixed_log_lh_only()`.
