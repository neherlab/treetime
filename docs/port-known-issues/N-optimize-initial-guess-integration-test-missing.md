# Integration test missing for initial guess mode dispatch

Unit tests cover `should_run_initial_guess()` in isolation, but no test exercises `run_optimize()` with explicit `InitialGuessMode::Always` or `InitialGuessMode::Never` values. A wiring bug in the conditional dispatch at [packages/treetime/src/commands/optimize/run.rs#L132-L136](../../packages/treetime/src/commands/optimize/run.rs#L132-L136) would go undetected.

## Current coverage

The golden master test `test_gm_optimize` exercises `run_optimize()` with default `Auto` mode. Since the test tree has all branch lengths, this exercises the "skip initial guess" path. An inverted condition would cause the golden master to diverge.

## Impact

Negligible. The unit tests verify the decision function. The golden master provides indirect integration coverage. The wiring is three lines of code.

## Proposed fix

Add integration tests constructing `TreetimeOptimizeArgs` with `InitialGuessMode::Never` and `InitialGuessMode::Always`, verifying that branch lengths are preserved vs overwritten respectively.
