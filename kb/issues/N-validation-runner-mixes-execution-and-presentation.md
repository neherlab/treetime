# Validation runner mixes execution, aggregation, and presentation

`TestRunner` defines six presentation methods that only forward to `ValidationConsole`; none of its three implementations overrides them [`packages/treetime-validation/src/testing/runners/runner.rs#L157-L213`](../../packages/treetime-validation/src/testing/runners/runner.rs#L157-L213).

`fn compute_all_metrics()` aggregates numerical results, applies unit conversions, and immediately formats every value as `String` [`packages/treetime-validation/src/testing/console/console_metrics.rs#L15`](../../packages/treetime-validation/src/testing/console/console_metrics.rs#L15). The runner orchestration also owns console output and plot, TSV, JSON, and summary publication.

## Consequences

- Runner implementations inherit presentation methods they neither customize nor need.
- Numeric metrics lose their types before JSON, TSV, console, or tests can consume them.
- Sorting or comparing rendered values requires parsing presentation strings.
- Adding an output format changes the same orchestration that owns test selection and algorithm execution.

## Required boundary

Keep `TestRunner` focused on suite/case selection and producing typed outcomes. Compute a typed aggregate metrics value with explicit units. Compose console, plot, TSV, JSON, and summary reporters around that value, formatting only at each renderer boundary.

Remove virtual presentation hooks unless a real runner implementation requires distinct presentation behavior.

## Validation

- All three runners implement only execution-related behavior.
- Typed aggregate tests compare numeric values and units without parsing strings.
- Renderer tests independently assert console, TSV, JSON, plot, and summary formatting.
- Existing validation artifacts retain their values and documented formatting.

## Proposed implementation plan

1. Remove the six `print_*` default methods of `trait TestRunner` (`runner.rs#L174-L213`); `fn run_tests_generic()` and `fn execute_single_test()` call `ValidationConsole` directly
2. Change `fn compute_all_metrics()` (`packages/treetime-validation/src/testing/console/console_metrics.rs`) to return a typed struct of aggregate metrics (numbers with their unit), and move the string formatting into the console, TSV and JSON writers
3. Move console output, plots, TSV, JSON and summary publication out of `run_tests_generic()` into reporters that take the typed outcomes and the aggregate
4. Tests: an aggregate test on hand-made `TestResult`s that compares the typed struct as a whole (oracle: hand calculation of means and unit conversions); one renderer test per format on a fixed aggregate. Validate with `just test-rs -E 'package(treetime-validation)'` and a run of each runner whose TSV and JSON outputs equal those before the change
