# Code quality and convention violations

## Summary

Convention violations across the codebase: fully qualified paths, inconsistent derives, manual loops replaceable by ndarray operations, module ordering violations, naming conventions, and miscellaneous quality issues.

## Instances

### ndarray over manual loops (7 instances)

1. Triple-nested loop replaceable by `sum_axis()` [packages/treetime/src/gtr/infer_gtr/common.rs#L100-L129](../../packages/treetime/src/gtr/infer_gtr/common.rs#L100-L129)
2. Manual outer product [packages/treetime/src/gtr/infer_gtr/common.rs#L59-L79](../../packages/treetime/src/gtr/infer_gtr/common.rs#L59-L79)
3. Manual W computation [packages/treetime/src/gtr/get_gtr.rs#L398-L409](../../packages/treetime/src/gtr/get_gtr.rs#L398-L409)
4. Manual einsum [packages/treetime/src/gtr/infer_gtr/site_specific.rs#L209-L223](../../packages/treetime/src/gtr/infer_gtr/site_specific.rs#L209-L223)
5. Zeros plus loop instead of `from_shape_fn` [packages/treetime-grid/src/interp_nonuniform.rs#L38-L53](../../packages/treetime-grid/src/interp_nonuniform.rs#L38-L53)
6. Manual swap loop instead of `reverse_inplace()` [packages/treetime-grid/src/grid_fn.rs#L397-L401](../../packages/treetime-grid/src/grid_fn.rs#L397-L401)
7. ndarray-to-`Vec` conversion [packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs#L43-L44](../../packages/treetime-validation/src/testing/metrics/aggregate/domain_agreement/error_stats.rs#L43-L44)

### NodeTimetree.nwk_comments() truncates date to 2 decimals

`fn NodeTimetree.nwk_comments()` [packages/treetime/src/payload/timetree.rs#L130-L136](../../packages/treetime/src/payload/timetree.rs#L130-L136)

Two decimal places gives ~3.6 days resolution, insufficient for fast-evolving pathogens sampled days apart.

### BrentOpt absolute epsilon on calendar time

Polytomy merge-time bracket [packages/treetime/src/timetree/optimization/polytomy.rs#L265](../../packages/treetime/src/timetree/optimization/polytomy.rs#L265)

`BrentOpt::new(parent_time + 1e-10, child_min_time - 1e-10)` uses absolute epsilon. For calendar times (e.g., 2020.5), 1e-10 is below f64 precision at that magnitude.

### edge_effective_length saturating_sub clamps to 0

- Sparse `fn edge_effective_length()` [packages/treetime/src/partition/marginal_sparse.rs#L255](../../packages/treetime/src/partition/marginal_sparse.rs#L255)
- Dense `fn edge_effective_length()` [packages/treetime/src/partition/marginal_dense.rs#L153](../../packages/treetime/src/partition/marginal_dense.rs#L153)

Returns 0 when non-char positions exceed sequence length. Used as denominator in `edge_subs().len() / edge_effective_length()`, producing division by zero.

### Brent cost function recomputes exp(eigvals\*t) per call

`mod method_brent` [packages/treetime/src/optimize/method_brent.rs](../../packages/treetime/src/optimize/method_brent.rs)

When multiple `OptimizationContribution` entries share the same GTR eigvals at the same branch length, the `exp_ev` vector is redundantly computed per contribution.

### Distribution::Empty fallback prevents polytomy merges

Polytomy distribution evaluation [packages/treetime/src/timetree/optimization/polytomy.rs#L216](../../packages/treetime/src/timetree/optimization/polytomy.rs#L216)

Silently prevents merges when branch-length distributions are Empty. Should log or error.

### TODO comments

1. Delimiter behavior [packages/treetime-io/src/concat.rs#L38-L42](../../packages/treetime-io/src/concat.rs#L38-L42)
2. Inefficient grid method [packages/treetime-grid/src/grid_fn.rs#L180](../../packages/treetime-grid/src/grid_fn.rs#L180)

### edge_effective_length clones gap vectors

- Dense `fn range_union()` call [packages/treetime/src/partition/marginal_dense.rs#L122](../../packages/treetime/src/partition/marginal_dense.rs#L122)
- Sparse `fn range_union()` call [packages/treetime/src/partition/marginal_sparse.rs#L380](../../packages/treetime/src/partition/marginal_sparse.rs#L380)

`fn range_union` clones both gap vectors on every call. Allocation pressure on heavily gapped alignments.

## Impact

- Division by zero in edge_effective_length
- Duplicate normalization code creates maintenance divergence risk

### Output conversion imports and enum formatting

Six server argument conversions place imports inside functions, obscuring their dependency surface:

- `impl From<ServerAncestralArgs> for TreetimeAncestralArgs` [packages/app-server/src/args.rs#L48-L54](../../packages/app-server/src/args.rs#L48-L54)
- `impl From<ServerClockArgs> for TreetimeClockArgs` [packages/app-server/src/args.rs#L132-L137](../../packages/app-server/src/args.rs#L132-L137)
- `impl From<ServerTimetreeArgs> for TreetimeTimetreeArgs` [packages/app-server/src/args.rs#L242-L249](../../packages/app-server/src/args.rs#L242-L249)
- `impl From<ServerMugrationArgs> for TreetimeMugrationArgs` [packages/app-server/src/args.rs#L358-L360](../../packages/app-server/src/args.rs#L358-L360)
- `impl From<ServerOptimizeArgs> for TreetimeOptimizeArgs` [packages/app-server/src/args.rs#L424-L430](../../packages/app-server/src/args.rs#L424-L430)
- `impl From<ServerPruneArgs> for TreetimePruneArgs` [packages/app-server/src/args.rs#L486-L490](../../packages/app-server/src/args.rs#L486-L490)

Two enums hand-write `Display`, duplicating their kebab-case CLI names:

- `impl Display for OutputSelection` [packages/treetime/src/commands/shared/output.rs#L163-L167](../../packages/treetime/src/commands/shared/output.rs#L163-L167)
- `impl Display for TopologyOrderTargetSourceArg` [packages/treetime/src/commands/shared/output.rs#L940-L947](../../packages/treetime/src/commands/shared/output.rs#L940-L947)

#### Axis A1: server conversion structure

- O1. Move the shared argument-type imports to module scope and retain the six explicit conversions. This exposes dependencies while preserving the current direct field mapping.
- O2. Extract focused conversion helpers for repeated shared argument groups, with their dependencies imported at module scope. This can also enforce identical defaults, but only where the conversion contracts are genuinely identical.

**Recommendation:** O1 for imports, then use O2 only for field groups proven identical by whole-structure tests. Avoid coupling command-specific defaults through a broad helper.

#### Axis A2: enum rendering

- O1. Derive `Display` with explicit kebab-case configuration from the same enum variants used by CLI parsing.
- O2. Retain manual implementations and add exhaustive tests tying every rendered value to the corresponding CLI name. This preserves bespoke control but duplicates the mapping.

**Recommendation:** O1. One derive-based mapping keeps CLI parsing and rendering exhaustive when variants are added.

## Related tickets

- [kb/tickets/convention-replace-manual-loops-with-ndarray-operations.md](../tickets/convention-replace-manual-loops-with-ndarray-operations.md)
- [kb/tickets/optimize-brent-cost-function-recomputes-exp-per-call.md](../tickets/optimize-brent-cost-function-recomputes-exp-per-call.md)
- [kb/tickets/optimize-edge-effective-length-saturating-sub-and-gap-clones.md](../tickets/optimize-edge-effective-length-saturating-sub-and-gap-clones.md)
