# GTR::new() panics on invalid dimensions

`GTR::new()` at [packages/treetime/src/gtr/gtr.rs#L76](../../packages/treetime/src/gtr/gtr.rs#L76) uses `assert!` for precondition validation: `assert!(abs(W.diag().sum()) < 1e-10)`. If the input rate matrix `W` has a non-zero diagonal sum, the function panics instead of returning a typed error.

`eig_single_site` (called from `GTR::new`) also uses `assert!` for dimension checks on the input matrices.

## Impact

`GTR::new()` returns `Result<Self, Report>`, so there is no API reason for it to panic. All callers already handle the `Result`. The panic is reachable from any code path that constructs a GTR from user-provided or inferred parameters with numerical noise that violates the precondition.

## Affected code

- Panic site: [packages/treetime/src/gtr/gtr.rs#L76](../../packages/treetime/src/gtr/gtr.rs#L76) -- `assert!(abs(W.diag().sum()) < 1e-10)`
- Also: `eig_single_site` at [packages/treetime/src/gtr/gtr.rs](../../packages/treetime/src/gtr/gtr.rs) uses `assert!` for preconditions

## Fix

Replace `assert!` calls with `make_error!` returns, using the existing `Result` return type. No API change needed since `GTR::new()` already returns `Result`.
