# Mechanical lint and quality fixes across codebase

Assorted convention violations and mechanical quality issues. Each item is independent and low-risk.

## Instances

### Fully qualified paths instead of imports (2 locations)

1. `packages/treetime-utils/src/io/compression.rs#L131` and `#L190` -- `std::io::Result<usize>` in trait impl signatures
2. `packages/treetime-utils/src/lib.rs#L20` -- `tikv_jemallocator::Jemalloc` in global static

### Inconsistent derives: SmartDefault vs std Default

`packages/treetime/src/commands/timetree/args.rs#L20` (and L29, L39, L49) -- three enums use `#[derive(Default)]` + `#[default]` (std), while `struct TreetimeTimetreeArgs` uses `SmartDefault`. Should use one approach consistently.

### Module ordering violations (4 locations)

1. `packages/treetime/src/commands/optimize/run.rs#L52` -- private fn before pub fn
2. `packages/treetime/src/commands/optimize/optimize_unified.rs#L729` -- pub fn after private fn
3. `packages/treetime/src/commands/prune/run.rs#L34` -- private fn before pub fn
4. `packages/treetime-ops/src/multiplication.rs#L10` -- private fn before pub fn

### Float-formatted BTreeMap keys

`packages/treetime-validation/src/testing/metrics/pointwise/structural.rs#L122-L126` -- `format!("{val:.12}")` used as `BTreeMap` key for numeric lookup. String-formatted floats as map keys are fragile.

### Section separator comment

`packages/treetime/src/commands/timetree/run.rs#L305` -- `// --- Postprocessing ---` banner comment. Should split file instead.

### Same-type tuple return

`packages/treetime-analytical/src/gaussian.rs#L30` -- `fn gaussian_product_params()` returns `(f64, f64, f64)`. A named struct would prevent argument-order bugs.

### Unchecked integer casts (3 locations)

1. `packages/treetime-cli/src/convert/usher.rs#L17` -- `i32 as usize` (negative wraps to large positive)
2. `packages/treetime-cli/src/convert/usher.rs#L34` -- `i32 as usize`
3. `packages/treetime-cli/src/convert/usher.rs#L56` -- `usize as i32` (truncation on large values)

### write_clock_regression_chart_svg bypasses file-helper conventions

`packages/treetime/src/cli/rtt_chart.rs#L27` -- passes `filepath` directly to `SVGBackend::new()` instead of using `fn create_file_or_stdout()` from `treetime-utils`.

### Normalize helpers duplicated

`packages/treetime/src/representation/partition/discrete.rs#L223` and `#L232` -- `fn normalize_inplace_1d` and `fn normalize_from_log_1d` are Array1 analogs of the Array2 versions in `marginal_dense.rs#L493` and `#L514`. The 1D versions add error handling that the 2D versions lack. Should be unified.

### #[allow(trivial_casts)] + magic damping 0.75

`packages/treetime/src/commands/timetree/run.rs#L477-L481` -- `let damping = 0.75;` is a hardcoded magic number (v0 default from `treeanc.py:1298`). No named constant or configuration.

### #[allow(clippy::useless_let_if_seq)]

`packages/treetime/src/commands/timetree/refinement.rs#L20` -- suppresses clippy lint. Could be rewritten as a single `let ... = if ...` expression.

### TODO/FIXME comments (6 locations)

1. `packages/treetime-io/src/nwk.rs#L73` -- TODO: parse nwk comments
2. `packages/treetime-io/src/concat.rs#L38-L42` -- 2 TODOs about delimiter behavior
3. `packages/treetime-cli/src/convert/auspice.rs#L73` and `#L80` -- 2 FIXMEs about assumed-zero values
4. `packages/treetime/src/commands/clock/run.rs#L64` -- TODO: empirical overdispersion value
5. `packages/treetime-grid/src/grid_fn.rs#L180` -- TODO: remove inefficient method

### PartitionFitchConfig::new() infallible Result

`packages/treetime/src/representation/partition/fitch_config.rs#L27` -- constructor wraps result in `Ok()` without validation. `Result` return type is misleading when the function cannot fail.

### calculate_diff_stats divides by zero on disjoint sets

`packages/treetime/examples/timetree_validation.rs#L411-L412` -- divides by `n` (count of overlapping keys). Disjoint node sets produce `n = 0`, yielding NaN.

## Related issues

- Source: [N-code-quality-conventions.md](../issues/N-code-quality-conventions.md) -- delete after full resolution
