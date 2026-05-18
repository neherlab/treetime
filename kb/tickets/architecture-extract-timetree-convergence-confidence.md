# Extract timetree convergence and confidence domain logic from commands

## Description

`commands/timetree/convergence/` and `commands/timetree/output/confidence.rs` contain domain logic that belongs in `src/timetree/`.

### Convergence domain logic

`commands/timetree/convergence/likelihood.rs` (94 lines): `compute_sequence_likelihood`, `compute_positional_likelihood`, `compute_coalescent_likelihood`. Pure domain computations, no CLI types.

`commands/timetree/convergence/sequence_changes.rs` (77 lines): `capture_ancestral_states`, `count_sequence_changes`, `AncestralStateSnapshot`. Pure domain computations, no CLI types.

`commands/timetree/convergence/metrics.rs`: `ConvergenceMetrics` struct (domain data) mixed with `TimetreeOptimizer` (CLI loop control) and `TreetimeOptimizerTraceCsvWriter` (I/O).

Move to `src/timetree/convergence/`:

- `likelihood.rs` (entire file)
- `sequence_changes.rs` (entire file)
- `ConvergenceMetrics` struct and `has_converged()` method

Keep in `commands/timetree/convergence/`:

- `TimetreeOptimizer` (loop control, consumes domain metrics)
- `TreetimeOptimizerTraceCsvWriter` (I/O)

### Confidence/rate-susceptibility domain logic

`commands/timetree/output/confidence.rs` (375 lines): rate susceptibility analysis, probit function, quadrature CI combination, CI extraction. ~60% domain logic (statistical computation), ~40% I/O and graph traversal.

Move to `src/timetree/confidence.rs`:

- `compute_rate_susceptibility`, `date_uncertainty_due_to_rate`, `combine_confidence`, `quantile_to_zscore`, `determine_rate_std`, `extract_confidence_intervals`, `NodeConfidenceInterval`
- Private helpers: `save_gammas`, `scale_gammas`, `collect_node_times`
- Constants: `CI_FRACTION`, `CI_LOWER_QUANTILE`, `CI_UPPER_QUANTILE`

Keep in `commands/timetree/output/confidence.rs`:

- `write_confidence_intervals` (TSV I/O)

## Validation

- `src/timetree/convergence/` compiles independently with no `commands/` imports
- `src/timetree/confidence.rs` compiles independently with no `commands/` imports
- All existing tests pass
- `commands/timetree/` imports from `crate::timetree::convergence` and `crate::timetree::confidence`

## Related issues

Source: [H-core-multi-client-architecture-library-purity](../issues/H-core-multi-client-architecture-library-purity.md)
