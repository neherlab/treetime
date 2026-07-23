# Validation runner mixes execution, aggregation, and presentation

`TestRunner` defines six presentation methods that only forward to `ValidationConsole`; none of its three implementations overrides them [`packages/treetime-validation/src/testing/runners/runner.rs#L19-L75`](../../packages/treetime-validation/src/testing/runners/runner.rs#L19-L75).

`fn compute_all_metrics()` aggregates numerical results, applies unit conversions, and immediately formats every value as `String` [`packages/treetime-validation/src/testing/console/console_metrics.rs#L16`](../../packages/treetime-validation/src/testing/console/console_metrics.rs#L16). The runner orchestration also owns console output and plot, TSV, JSON, and summary publication.

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
