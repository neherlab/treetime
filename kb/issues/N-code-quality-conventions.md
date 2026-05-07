# Code quality and convention violations

## Summary

Convention violations across the codebase: fully qualified paths, inconsistent derives, manual loops replaceable by ndarray operations, module ordering violations, naming conventions, and miscellaneous quality issues.

## Instances

### Fully qualified paths instead of imports (2 locations)

1. `treetime-utils/src/io/compression.rs:131,190:` `std::io::Result<usize>` in trait impl signatures
2. `treetime-utils/src/lib.rs:20:` `tikv_jemallocator::Jemalloc` in global static

### Inconsistent derives: SmartDefault vs std Default

`packages/treetime/src/commands/timetree/args.rs:20,29,39,49:`

Three enums use `#[derive(Default)]` + `#[default]` (std), while `struct TreetimeTimetreeArgs` uses `SmartDefault`. Should use one approach consistently.

### ndarray over manual loops (7 instances)

1. `gtr/infer_gtr/dense.rs:100-129:` triple-nested loop replaceable by `sum_axis()`
2. `gtr/infer_gtr/dense.rs:59-79:` manual outer product
3. `gtr/get_gtr.rs:398-409:` manual W computation
4. `gtr/infer_gtr/site_specific.rs:209-223:` manual einsum
5. `treetime-grid/src/interp_nonuniform.rs:38-53:` zeros + loop instead of `from_shape_fn`
6. `treetime-grid/src/grid_fn.rs:397-401:` manual swap loop instead of `reverse_inplace()`
7. `treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs:43-44:` ndarray to Vec conversion

### Module ordering violations (4 instances)

1. `commands/optimize/run.rs:52:` private fn before pub fn
2. `commands/optimize/optimize_unified.rs:729:` pub fn after private fn
3. `commands/prune/run.rs:34:` private fn before pub fn
4. `treetime-ops/src/multiplication.rs:10:` private fn before pub fn

### Float-formatted BTreeMap keys

`packages/treetime-validation/src/testing/metrics/pointwise/structural.rs:122-126:`

`format!("{val:.12}")` used as `BTreeMap` key for numeric lookup. String-formatted floats as map keys are fragile (formatting changes break lookups).

### Section separator comment

`packages/treetime/src/commands/timetree/run.rs:305:` `// --- Postprocessing ---` banner comment. Should split file instead.

### Same-type tuple return

`packages/treetime-analytical/src/gaussian.rs:30:`

`fn gaussian_product_params()` returns `(f64, f64, f64)`. A named struct would prevent argument-order bugs.

### Unchecked integer casts (3 instances)

- `treetime-cli/src/convert/usher.rs:17:` `i32 as usize` (negative wraps to large positive)
- `treetime-cli/src/convert/usher.rs:34:` `i32 as usize`
- `treetime-cli/src/convert/usher.rs:56:` `usize as i32` (truncation on large values)

### Alphabet serde skip fields not rebuilt

`packages/treetime/src/alphabet/alphabet.rs:48-56:`

Four `#[serde(skip)]` fields with no post-deserialization rebuild. If `struct Alphabet` is deserialized, lookup tables are empty/default.

### Missing empty-canonical rejection

`packages/treetime/src/alphabet/alphabet_config.rs:72:`

`fn validate()` does not reject an empty canonical character set. An empty alphabet produces division-by-zero in profile construction.

### calculate_diff_stats divides by zero on disjoint sets

`packages/treetime/examples/timetree_validation.rs:411-412:`

Divides by `n` (count of overlapping keys). Disjoint node sets produce `n = 0`, yielding NaN.

### write_clock_regression_chart_svg bypasses file-helper conventions

`packages/treetime/src/cli/rtt_chart.rs:27:`

Passes `filepath` directly to `SVGBackend::new()` instead of using `fn create_file_or_stdout()` from `treetime-utils`.

### BrentOpt absolute epsilon on calendar time

`packages/treetime/src/commands/timetree/optimization/polytomy.rs:265:`

`BrentOpt::new(parent_time + 1e-10, child_min_time - 1e-10)` uses absolute epsilon. For calendar times (e.g., 2020.5), 1e-10 is below f64 precision at that magnitude.

### Normalize helpers duplicated in discrete.rs

`packages/treetime/src/representation/partition/discrete.rs:223,232:`

`fn normalize_inplace_1d` and `fn normalize_from_log_1d` are Array1 analogs of the Array2 versions in `marginal_dense.rs:493,514:`. The 1D versions add error handling that the 2D versions lack. Should be unified.

### edge_effective_length saturating_sub clamps to 0

`packages/treetime/src/representation/partition/marginal_sparse.rs:255:` and `marginal_dense.rs:153:`

Returns 0 when non-char positions exceed sequence length. Used as denominator in `edge_subs().len() / edge_effective_length()`, producing division by zero.

### Brent cost function recomputes exp(eigvals\*t) per call

`packages/treetime/src/commands/optimize/method_brent.rs`

When multiple `OptimizationContribution` entries share the same GTR eigvals at the same branch length, the `exp_ev` vector is redundantly computed per contribution.

### #[allow(trivial_casts)] + magic damping 0.75

`packages/treetime/src/commands/timetree/run.rs:477-481:`

`let damping = 0.75;` is a hardcoded magic number (v0 default from `treeanc.py:1298`). No named constant or configuration.

### #[allow(clippy::useless_let_if_seq)]

`packages/treetime/src/commands/timetree/refinement.rs:20:`

Suppresses clippy lint. Could be rewritten as a single `let ... = if ...` expression.

### Distribution::Empty fallback prevents polytomy merges

`packages/treetime/src/commands/timetree/optimization/polytomy.rs:216:`

Silently prevents merges when branch-length distributions are Empty. Should log or error.

### TODO/FIXME comments (6 locations)

1. `treetime-io/src/nwk.rs:73:` TODO: parse nwk comments
2. `treetime-io/src/concat.rs:38-42:` 2 TODOs about delimiter behavior
3. `treetime-cli/src/convert/auspice.rs:73,80:` 2 FIXMEs about assumed-zero values
4. `commands/clock/run.rs:64:` TODO: empirical overdispersion value
5. `treetime-grid/src/grid_fn.rs:180:` TODO: remove inefficient method

### edge_effective_length clones gap vectors

`packages/treetime/src/representation/partition/marginal_dense.rs:122:` and `marginal_sparse.rs:380:`

`fn range_union` clones both gap vectors on every call. Allocation pressure on heavily gapped alignments.

### PartitionFitchConfig::new() infallible Result

`packages/treetime/src/representation/partition/fitch_config.rs:27:`

Constructor wraps result in `Ok()` without validation. `Result` return type is misleading when the function cannot fail.

## Impact

- Integer cast bugs produce wrong indices on negative or large values
- Division by zero in edge_effective_length and diff_stats
- Magic numbers without named constants reduce maintainability
- Duplicate normalization code creates maintenance divergence risk
