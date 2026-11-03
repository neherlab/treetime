# Grid point limit tests assume a 64-bit `usize`

`MaxGridPoints::MAX` at [packages/treetime-grid/src/max_grid_points.rs#L24](../../packages/treetime-grid/src/max_grid_points.rs#L24) is the smaller of $2^{53} - 1$ and `usize::MAX`, so the production code is correct on 32-bit targets. The tests are not:

- [packages/treetime-grid/src/__tests__/grid.rs#L95-L96](../../packages/treetime-grid/src/__tests__/grid.rs#L95-L96): the cases `inexact_point_count` and `first_inexact_count` pass `(1 << 53) + 2` and `1 << 53` as `usize`
- [packages/treetime-grid/src/__tests__/max_grid_points.rs#L27](../../packages/treetime-grid/src/__tests__/max_grid_points.rs#L27): the case `maximum` expects `(1 << 53) - 1`

On a target with a 32-bit `usize`, such as `wasm32-unknown-emscripten`, these literals overflow and the test crate does not compile. The expected error texts and the schema bound also name the 64-bit value. No recipe compiles the tests for a 32-bit target today: `just lint-wasm` and `just lint-wasm-mt` check only the binaries of `app-wasm`, and `just test-wasm` runs the built module.

## Fix

Derive the expected values of these cases from `MaxGridPoints::MAX`, or gate the cases that need $2^{53}$ points on `target_pointer_width = "64"`, when a recipe starts to compile tests for a 32-bit target.
