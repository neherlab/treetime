# Small profile products spend most of their BLAS time on OpenBLAS packing

The dense passes multiply an L-by-s profile by an s-by-s matrix (s = 4 for nucleotides, 20 for amino acids). ndarray sends every matrix product with a dimension above 7 to `cblas_dgemm`, so these products run through the full OpenBLAS general matrix multiplication. OpenBLAS first copies ("packs") the profile into a work buffer, and for s = 4 this copy costs 2.5 to 3 times as much as the arithmetic.

## Measurement

Profiling build of `rust` `abad3947e`, on the build host (i9-10900K, OpenBLAS Haswell kernels, 20 threads), in a window with no other heavy jobs. Share of all CPU samples (`perf record`) and totals from uprobes on the same binary (bpftrace, times summed over threads):

| Run | Packing (`dgemm_*copy`) | Arithmetic (`dgemm_kernel`) | `cblas_dgemm` calls | Time in `cblas_dgemm` | Time in the buffer pool | Sleep on its lock |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `ancestral --representation=dense` `rsv/a/2000` | 7.43% | 2.65% | 7990 | 9.49 s | 72 ms | 8 ms |
| `optimize --representation=dense` `rsv/a/500` | 6.65% | 2.17% | 12350 | 9.40 s | 75 ms | 2.5 ms |
| `ancestral --representation=dense` `flu/h3n2/500` | 0.90% | 2.49% | 1832 | 0.09 s | 6 ms | 1.5 ms |

- About half of the time in `cblas_dgemm` is not BLAS work: it is the first-touch page faults of the new output arrays, which `dgemm_beta` zeroes ([N-ancestral-dense-messages-fill-memory-through-first-touch-page-faults.md](N-ancestral-dense-messages-fill-memory-through-first-touch-page-faults.md))
- The buffer pool of OpenBLAS (`blas_memory_alloc`, `blas_memory_free`, one global mutex) takes under 1 percent of the time in `cblas_dgemm`, and threads sleep on its lock for less than 10 ms per run. The threaded OpenBLAS build with per-thread buffers (`USE_TLS=1`) would remove only this cost
- Removing packing and the call overhead saves at most about 7 percent of CPU time on these dense runs

## Location

- `const GEMM_BLAS_CUTOFF: usize = 7` in ndarray 0.17.2 (`src/linalg/impl_linalg.rs`): a product with any dimension above 7 goes to `cblas_dgemm`
- The products: `fn GTR.evolve()` [packages/treetime/src/gtr/gtr.rs#L154](../../packages/treetime/src/gtr/gtr.rs#L154), `fn GTR.propagate_profile()` [packages/treetime/src/gtr/gtr.rs#L174](../../packages/treetime/src/gtr/gtr.rs#L174), `fn in_eigenbasis()` [packages/treetime/src/gtr/gtr.rs#L259-L267](../../packages/treetime/src/gtr/gtr.rs#L259-L267), and `fn get_coefficients()` [packages/treetime/src/partition/optimize/dense.rs#L5-L12](../../packages/treetime/src/partition/optimize/dense.rs#L5-L12)
- OpenBLAS (commit `de5d4914`) has a path without packing and without the buffer pool for small products (`M * N * K <= 1e6`), but only the SkylakeX-class AVX-512 kernels and the ARM SVE kernels enable it. The Haswell and Zen kernels use `kernel/generic/gemm_small_matrix_permit.c`, which always returns 0. The gain of a change therefore differs by CPU

## Fix direction

Compute products with a small s without OpenBLAS, so that no packing copy, buffer pool, or dispatch runs per call. Keep OpenBLAS for large products and for LAPACK.

> [!IMPORTANT]
> **Decision required.** Each option conflicts with a project rule or adds maintenance:
>
> - Call `matrixmultiply::dgemm` (the pure-Rust backend of ndarray, already in the dependency tree) at these call sites. The "Tech stack" rule names only ndarray and ndarray-linalg on OpenBLAS for linear algebra, so it needs an entry. matrixmultiply also packs, but allocates its buffer through jemalloc's per-thread cache and runs no global lock
> - Write a kernel for a fixed small s with ndarray views. The rules rank custom helpers last
> - Change the dispatch rule of ndarray upstream, for example to use BLAS only when the smallest dimension is above the cutoff
>
> Measure the chosen option before and after on the three runs above, at 1 and 20 threads, on the same host.
