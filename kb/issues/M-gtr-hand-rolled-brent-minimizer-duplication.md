# GTR rate optimization uses hand-rolled Brent instead of argmin wrapper

`gtr/refinement.rs` contains a 100-line manual Brent minimizer (`brent_minimize_bracketed`) for GTR rate optimization. A production-quality argmin-based Brent wrapper already exists in `optimize/method_brent.rs` with proper error handling via `Result<f64, Report>`, machine-epsilon tolerance, and sqrt-space transform support (`brent_sqrt_inner`).

The hand-rolled version returns a bare `f64` with no error propagation, uses a fixed tolerance of `1e-8` vs machine epsilon, allows 100 iterations vs 50, and silently returns the best guess on exhaustion.

## Locations

- `packages/treetime/src/gtr/refinement.rs:147-247:` `fn brent_minimize_bracketed()` -- hand-rolled 100-line Brent with golden-section fallback
- `packages/treetime/src/gtr/refinement.rs:130:` sole production caller, in `optimize_gtr_rate()`, operating in sqrt(mu) space
- `packages/treetime/src/optimize/method_brent.rs:41-67:` `fn brent_inner()` -- argmin `BrentOpt` wrapper with `Result<f64, Report>`
- `packages/treetime/src/optimize/method_brent.rs:85-114:` `fn brent_sqrt_inner()` -- same wrapper with sqrt-space transform (mirrors `optimize_gtr_rate` usage)

## Impact

Two independent implementations of the same algorithm with divergent error handling and tolerance. Fixes to one do not propagate to the other. The hand-rolled version's silent failure mode masks optimization failures.

## Action

Replace `brent_minimize_bracketed` with the argmin wrapper. The sqrt-space transform in `optimize_gtr_rate` maps directly to `brent_sqrt_inner`. Requires adapting the cost function interface from `FnMut(f64) -> f64` to argmin's `CostFunction` trait.
