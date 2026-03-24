# Optimize loop not extracted from I/O wrapper

`run_optimize()` bundles file I/O (read tree, read alignment, write output) with the core optimization loop (setup partitions, iterate marginal reconstruction + branch length optimization + damping). Tests that exercise the optimization loop must reimplement it rather than calling the production code path.

## Problem

Tests in `test_gm_optimize.rs`, `test_damping.rs`, and `test_convergence/test_convergence_iterations.rs` each contain their own copy of the optimization loop. If `run_optimize()` changes (new steps, different ordering, additional setup), these test copies diverge silently from production.

## Affected tests

- `packages/treetime/src/commands/optimize/__tests__/test_gm_optimize.rs` - `setup_and_run()` helper reimplements the full loop
- `packages/treetime/src/commands/optimize/__tests__/test_damping.rs` - `test_damped_optimization_converges()`, `test_damped_optimization_does_not_regress()`
- `packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_iterations.rs` - all three tests
- `packages/treetime/src/commands/optimize/__tests__/test_convergence/test_convergence_support.rs` - `setup_partitions()` duplicates partition creation

## Solution

Extract the core optimization loop from `run_optimize()` into a testable function that takes in-memory data structures and returns results. `run_optimize()` becomes an I/O wrapper that calls the extracted function. Tests call the extracted function directly.
