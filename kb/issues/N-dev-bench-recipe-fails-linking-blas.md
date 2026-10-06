# `just bench` fails to link the test harness of treetime-utils

`just bench -- <filter>` runs `cargo bench --locked --workspace --benches` with the shipped dist flags ([justfile](../../justfile), recipe `bench`). `--benches` also builds every library whose `bench` flag is set, which is the default for libraries, so cargo builds the test harness of `treetime-utils` in the bench profile. Linking it fails in the development container:

```text
rust-lld: error: undefined symbol: cblas_dgemm
rust-lld: error: undefined symbol: cblas_dgemv
rust-lld: error: undefined symbol: cblas_ddot
error: could not compile `treetime-utils` (lib test) due to 1 previous error
```

The symbols come from ndarray's `dot()`, which uses BLAS, in tests such as `test_matvec_3d_matches_per_site_dot` in [packages/treetime-utils/src/array/__tests__/test_batched.rs](../../packages/treetime-utils/src/array/__tests__/test_batched.rs). The dev and release test builds link them.

A single benchmark builds and runs without the recipe: `./dev/docker/run just build --profile bench -p <crate> --bench <name>`, then the built binary with `--bench`.

> [!IMPORTANT]
> **Investigation required.** Reproduce the failure on the tip of `rust` and find why the bench profile with the dist flags does not link OpenBLAS into library test harnesses. Possible fixes: set `bench = false` on the libraries without benchmarks, or link OpenBLAS in that profile as the test profiles do.
