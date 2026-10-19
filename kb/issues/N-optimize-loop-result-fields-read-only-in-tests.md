# Two fields of the optimization loop result are read only in tests

`struct OptimizeLoopResult` keeps the log-likelihood of every iteration (`lh_history`) and the iteration and reason of the stop (`stopped_at`), but only test code reads them; production builds carry `#[cfg_attr(not(test), expect(dead_code, reason = "read only in test code"))]` on both fields [packages/treetime/src/optimize/run_loop.rs#L115-L131](../../packages/treetime/src/optimize/run_loop.rs#L115-L131). The loop computes both for its own stop rule [packages/treetime/src/optimize/run_loop.rs#L49-L68](../../packages/treetime/src/optimize/run_loop.rs#L49-L68), so storing them costs little, but a result field that no command, output or app reads is a second channel of the run state beside `ConvergenceSummary`, and its tests check data that users never see.

## Fix direction

Either report the stop reason and the log-likelihood trace to the user (the run log or the trace output of `treetime optimize`), so the fields have a production reader, or keep them out of the result and let the tests observe the loop through its public outputs.
