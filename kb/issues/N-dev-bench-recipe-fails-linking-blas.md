# `just bench` fails to link the test harness of treetime-utils

`just bench -- <filter>` runs `cargo bench --locked --workspace --benches` with the shipped dist flags ([justfile#L609-L610](../../justfile#L609-L610), recipe `bench`). `--benches` also builds every library whose `bench` flag is set, which is the default for libraries, so cargo builds the test harness of `treetime-utils` in the bench profile. Linking it fails in the development container:

```text
rust-lld: error: undefined symbol: cblas_dgemm
rust-lld: error: undefined symbol: cblas_dgemv
rust-lld: error: undefined symbol: cblas_ddot
error: could not compile `treetime-utils` (lib test) due to 1 previous error
```

## Cause

`export_dist_flags()` in [dev/lib/dist-flags.sh](../../dev/lib/dist-flags.sh) sets `RUSTFLAGS` to the CPU flags only. `RUSTFLAGS` replaces the target `rustflags` of [.cargo/config.toml#L34-L41](../../.cargo/config.toml#L34-L41), which add `-l:libopenblas.a` and the Fortran runtime to every link line. rustc links the `openblas-src` library only into targets that reach a crate naming it (`ndarray-linalg`). The test harness of `treetime-utils` calls BLAS only through ndarray: `outer()` [packages/treetime-utils/src/array/ndarray.rs#L20-L24](../../packages/treetime-utils/src/array/ndarray.rs#L20-L24) calls `dot()`, which ndarray sends to `cblas_dgemm`. With the dist flags, its link line has no OpenBLAS. The dev and release test builds keep the target `rustflags` and link it.

Only `treetime` and `util-newick` declare benchmarks (`[[bench]]` in their `Cargo.toml`).

A single benchmark builds and runs without the recipe: `./dev/docker/run just build --profile bench -p <crate> --bench <name>`, then the built binary with `--bench`.

## Possible fixes

- Add the OpenBLAS and Fortran runtime link arguments to `RUSTFLAGS` in `export_dist_flags()` for the native target, so every dist, profiling, and bench build links them as the test profiles do
- Set `bench = false` on the libraries without benchmarks, so `--benches` does not build their test harnesses
