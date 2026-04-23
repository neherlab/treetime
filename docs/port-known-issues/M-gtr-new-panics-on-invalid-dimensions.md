# GTR::new() panics on invalid dimensions

`GTR::new()` at [packages/treetime/src/gtr/gtr.rs#L231](../../packages/treetime/src/gtr/gtr.rs#L231) uses `assert_eq!` for dimension validation at L234-L246. `eig_single_site` (called from `GTR::new`) at [packages/treetime/src/gtr/gtr.rs#L76](../../packages/treetime/src/gtr/gtr.rs#L76) uses `assert!(abs(W.diag().sum()) < 1e-10)` for diagonal sum validation. If the input rate matrix `W` has a non-zero diagonal sum, the function panics instead of returning a typed error.

## Impact

`GTR::new()` returns `Result<Self, Report>`, so there is no API reason for it to panic. All callers already handle the `Result`. The panic is reachable from any code path that constructs a GTR from user-provided or inferred parameters with numerical noise that violates the precondition.

## Affected code

- Dimension asserts: [packages/treetime/src/gtr/gtr.rs#L234-L246](../../packages/treetime/src/gtr/gtr.rs#L234-L246) -- `assert_eq!` in `GTR::new()`
- Diagonal sum assert: [packages/treetime/src/gtr/gtr.rs#L76](../../packages/treetime/src/gtr/gtr.rs#L76) -- `assert!(abs(W.diag().sum()) < 1e-10)` in `eig_single_site`

## Fix

Replace `assert!` calls with `make_error!` returns, using the existing `Result` return type. No API change needed since `GTR::new()` already returns `Result`.
