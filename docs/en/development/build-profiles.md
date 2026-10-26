---
title: Build profiles
description: The cargo profiles of development, release, and profiling builds.
---

- `dev` (`just build`, `just run`, the tests): unoptimized workspace crates, dependencies at `opt-level = 2`, full debug info. Every build for `x86_64-unknown-linux-gnu`, the development builds included, compiles Rust and C code for the CPU level of the shipped x86_64 binaries, x86-64-v2 (`.cargo/config.toml`)
- `dev-opt` (`just run-dev-opt`): `dev` with optimized workspace crates, for long runs on real datasets; rebuilds are slower
- `release` (`just build-release`, `just run-release`, `just example`, `just smoke`): optimized and fast to rebuild, without LTO
- `dist` (`just build-dist`, `just cross`, `just build-napi dist`, `just cross-napi`, nightly releases): the shipped binary and Node addon, with fat LTO, one codegen unit, and the CPU and link flags of the shipped builds from `export_dist_flags` of `dev/lib/dist-flags.sh` (`-C target-cpu=x86-64-v2` on x86-64, no frame pointers). These flags replace the target rustflags of the Cargo config files, so the native dist builds link with rust-lld. The glibc binaries of `build-dist`, `build-profiling`, `build-profiling-heap`, and `bench` link the GCC unwinder statically and need no `libgcc_s.so.1`, as the cross builds do; the Node addon keeps the shared unwinder of the Electron process
- `profiling` (`just build-profiling`, `just profile`): `dist` with full debug info; without frame pointers, profilers unwind call stacks with DWARF (`dev/profile` does) or LBR
- `profiling-heap` (`just build-profiling-heap`): `profiling` with jemalloc heap profiling compiled in, for heap profiles only
- `bench` (`just bench`): the `dist` settings
