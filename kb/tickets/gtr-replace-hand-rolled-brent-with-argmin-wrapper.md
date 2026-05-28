# Replace hand-rolled Brent minimizer with argmin wrapper in GTR rate optimization

`optimize_gtr_rate()` uses a 100-line manual Brent implementation with no error handling. The existing argmin-based `brent_sqrt_inner` already supports the same sqrt-space transform.

## Current state

`gtr/refinement.rs:147-247` contains `brent_minimize_bracketed()` returning bare `f64`. Called at line 130 in sqrt(mu) space with tolerance `1e-8` and 100 iterations.

## Target state

`optimize_gtr_rate()` uses `brent_sqrt_inner` (or a new variant) from `optimize/method_brent.rs`, gaining `Result<f64, Report>` error propagation and machine-epsilon tolerance. `brent_minimize_bracketed` is deleted.

## Implementation

1. Adapt `optimize_gtr_rate`'s cost function to implement argmin's `CostFunction` trait (or use the existing generic interface in `method_brent.rs`)
2. Replace the `brent_minimize_bracketed` call at `refinement.rs:130` with `brent_sqrt_inner` or equivalent
3. Propagate the `Result` return type through `optimize_gtr_rate` callers
4. Delete `brent_minimize_bracketed` and its tests (move any valuable test cases to `optimize/__tests__/`)
5. Verify GTR rate optimization produces identical results on test datasets

## Related issues

Source: `kb/issues/M-gtr-hand-rolled-brent-minimizer-duplication.md`
