---
title: Performance
description: Benchmarks, profiles, the allocator, and heap profiles.
---

- `just bench` runs the benchmarks
- `just build-dist treetime` builds the shipped binary into `.out/treetime`; use it with `hyperfine`, which the container provides
- `just profile treetime -- <args>` (host) samples a profile with samply or perf; read `dev/profile --help` first

## Allocator

On x86_64 and aarch64 Linux, the CLI and the API server allocate through jemalloc. jemalloc is built without a symbol prefix (the `override_allocator_on_supported_platforms` feature of `tikv-jemalloc-sys`), so the binary defines `malloc`, `free`, and the other allocator functions of the C library itself. Rust code, the C libraries (OpenBLAS, liblzma, libzstd), and the C library all allocate through jemalloc. The glibc binaries export these functions, so glibc binds its own calls to them. The static musl binaries keep the musl allocator only for the internal allocations that musl also frees itself. Other targets and the Node addon use the allocator of the system.

- `MALLOC_CONF` configures jemalloc, for example `MALLOC_CONF=dirty_decay_ms:1000`. The settings of the glibc allocator (`MALLOC_ARENA_MAX`, the `glibc.malloc` tunables) have no effect
- `dev/cross/check` fails when a Linux binary does not take these functions from jemalloc

## Heap profiles

Heap profilers that replace `malloc` through `LD_PRELOAD` (heaptrack) do not see the allocations of the CLI, because the binary binds its calls to its own `malloc`. jemalloc's own heap profiler sees them. It is compiled into the `profiling-heap` build only, because it changes the allocator code: take CPU profiles and timings from the `profiling` and `dist` builds.

```bash
./dev/docker/run just build-profiling-heap
mkdir -p tmp/heap
MALLOC_CONF=prof:true,prof_gdump:true,prof_final:true,lg_prof_sample:17,prof_prefix:tmp/heap/treetime \
  .build/container/cargo/profiling-heap/treetime <command> <args>
peak="$(ls -v tmp/heap/*.u*.heap | tail -n 1)"
jeprof --collapsed .build/container/cargo/profiling-heap/treetime "${peak}" > tmp/heap/peak.folded
flamegraph.pl --countname=bytes tmp/heap/peak.folded > tmp/heap/peak.svg
```

- `prof_gdump` writes a `.u<n>.heap` dump each time the heap grows past its previous maximum, so the dump with the highest number is the closest to the peak. A run on a few hundred sequences writes thousands of dumps and several hundred MB; delete the others once the peak dump is chosen. `prof_final` writes one more dump (`.f.heap`) at exit
- `lg_prof_sample:17` samples once per 128 KiB allocated on average (the default is 512 KiB)
- `--collapsed` keeps full Rust paths and inlined frames. `--text` shortens Rust names to their last segment and charges every byte to jemalloc's sampler, which makes it hard to read
- `jeprof --base=<earlier>.heap` subtracts an earlier dump, and `--alloc_space` reports total allocation volume when the run had `prof_accum:true`
- `jeprof` comes with the `libjemalloc-dev` package on Debian and Ubuntu; `flamegraph.pl` is in [FlameGraph](https://github.com/brendangregg/FlameGraph)
