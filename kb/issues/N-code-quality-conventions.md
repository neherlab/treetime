# Code quality and convention violations

## Summary

Convention violations across the codebase: fully qualified paths, inconsistent derives, manual loops replaceable by ndarray operations, module ordering violations, naming conventions, and miscellaneous quality issues.

## Instances

### ndarray over manual loops (7 instances)

1. `gtr/infer_gtr/dense.rs:100-129:` triple-nested loop replaceable by `sum_axis()`
2. `gtr/infer_gtr/dense.rs:59-79:` manual outer product
3. `gtr/get_gtr.rs:398-409:` manual W computation
4. `gtr/infer_gtr/site_specific.rs:209-223:` manual einsum
5. `treetime-grid/src/interp_nonuniform.rs:38-53:` zeros + loop instead of `from_shape_fn`
6. `treetime-grid/src/grid_fn.rs:397-401:` manual swap loop instead of `reverse_inplace()`
7. `treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs:43-44:` ndarray to Vec conversion

### NodeTimetree.nwk_comments() truncates date to 2 decimals

`packages/treetime/src/payload/timetree.rs:135:`

Two decimal places gives ~3.6 days resolution, insufficient for fast-evolving pathogens sampled days apart.

### BrentOpt absolute epsilon on calendar time

`packages/treetime/src/commands/timetree/optimization/polytomy.rs:265:`

`BrentOpt::new(parent_time + 1e-10, child_min_time - 1e-10)` uses absolute epsilon. For calendar times (e.g., 2020.5), 1e-10 is below f64 precision at that magnitude.

### edge_effective_length saturating_sub clamps to 0

`packages/treetime/src/partition/marginal_sparse.rs:255:` and `marginal_dense.rs:153:`

Returns 0 when non-char positions exceed sequence length. Used as denominator in `edge_subs().len() / edge_effective_length()`, producing division by zero.

### Brent cost function recomputes exp(eigvals\*t) per call

`packages/treetime/src/commands/optimize/method_brent.rs`

When multiple `OptimizationContribution` entries share the same GTR eigvals at the same branch length, the `exp_ev` vector is redundantly computed per contribution.

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

`packages/treetime/src/partition/marginal_dense.rs:122:` and `marginal_sparse.rs:380:`

`fn range_union` clones both gap vectors on every call. Allocation pressure on heavily gapped alignments.

## Impact

- Division by zero in edge_effective_length
- Duplicate normalization code creates maintenance divergence risk
