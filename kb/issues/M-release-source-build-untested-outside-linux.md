# Source build untested outside Linux

## Summary

A source build outside the project containers links the OpenBLAS of the system through pkg-config (the build script of `packages/treetime-openblas`). This works on Linux. On macOS with Homebrew OpenBLAS, on FreeBSD, and on Windows with the OpenBLAS of MSYS2 or conda-forge, no source build has run. Packagers (bioconda, conda-forge, Homebrew, Nix) build from source with their own OpenBLAS, so the package channels for these systems depend on this path.

## Evidence

- `packages/treetime-openblas/build.rs` names the install commands for Debian and Ubuntu, macOS (Homebrew), FreeBSD, and Windows (MSYS2, conda-forge). Only the Linux commands were tested
- Linux, verified on Ubuntu 24.04 without mise, the project archives, or mold, with `build-essential`, `pkg-config`, `protobuf-compiler`, and the Rust toolchain of `rust-toolchain.toml`:
  - Without `libopenblas-dev`, `cargo build --locked -p app-cli` stops in the build script of `treetime-openblas` with the text of pkg-config, the install commands of each platform, and a pointer to the section "Host" of `docs/en/development/setup.md`
  - With `libopenblas-dev`, pkg-config finds `openblas.pc` without extra `PKG_CONFIG_PATH` (through the alternatives link of Ubuntu; `-L/usr/lib/x86_64-linux-gnu/openblas-pthread/ -lopenblas`). The build passes, the binary links `libopenblas.so.0` and `libgfortran.so.5`, `treetime debug` reports OpenBLAS 0.3.26 with 16 threads on a 16-CPU machine, and `treetime --version` and `ancestral` on `data/flu/h3n2/20` run
  - With mise 2026.9.10 after `mise install openblas liblzma`, `mise exec -- cargo build --locked -p app-cli` gives a binary that links no OpenBLAS and no liblzma (`ldd` lists `libgcc_s`, `libm`, `libc`), and `treetime debug` reports the single-threaded `DYNAMIC_ARCH` OpenBLAS 0.3.34 of the project archive
  - The UShER protobuf schemas (`packages/util-usher-mat/build.rs`, `prost-build`) need `protoc`, which mise provides and which a build without mise takes from the system

> [!IMPORTANT]
> **Investigation required.** Build the CLI from source on macOS with Homebrew OpenBLAS, on FreeBSD, and on Windows with MSYS2 and with conda-forge, and record what fails (pkg-config lookup, link names, LAPACK symbols, the Fortran runtime).

## Impact

- Package recipes for macOS, FreeBSD, and Windows cannot rely on the source build until it has run there

## Potential solutions

- Run the source build once on each of these systems, correct the install commands of `packages/treetime-openblas/build.rs` and the section "Host" of `docs/en/development/setup.md` where they fail, and record the result here
