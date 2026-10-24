# timetree uses 5 to 10 times more CPU time at 20 threads than at 1 thread

At 20 worker threads, `timetree` uses 5.5 to 9.5 times the CPU time of a run at 1 thread, while its wall time drops only by a factor of 2.0 to 2.9. Shared physical cores and a lower all-core clock frequency explain at most about 2.5 times of this. The remaining CPU time is spent in real work and not in idle or scheduler code, so part of the work is probably done more than once when the run uses more threads.

## Measurement

Host with 10 physical cores and 2 hardware threads per core (SMT on). User CPU time and wall time in seconds, median of 3 interleaved runs, before and after the tree scans replaced the graph pass:

| Case | Revision | CPU, 1 thread | CPU, 20 threads | CPU ratio | Wall, 1 thread | Wall, 20 threads | Wall speedup |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `timetree rsv/a/2000` | `602a051bb` (graph pass) | 46.6 | 258.3 | 5.5 | 49.6 | 18.6 | 2.7 |
| `timetree rsv/a/2000` | tree scans | 45.5 | 287.9 | 6.3 | 48.3 | 16.6 | 2.9 |
| `timetree dengue/2000` | `602a051bb` (graph pass) | 119.7 | 1021.8 | 8.5 | 129.0 | 64.9 | 2.0 |
| `timetree dengue/2000` | tree scans | 117.6 | 1119.2 | 9.5 | 126.6 | 61.8 | 2.0 |

- Both revisions have the low efficiency. At 1 thread the tree scans use 2 percent less CPU time; at 20 threads they keep more threads busy and use about 10 percent more
- The CPU profiles at 20 threads (`perf record`) show real-work symbols at the top: FFT planning with about 5 percent of the samples, interpolation, `exp`, and `sincos`. Idle and scheduler symbols have less than 0.3 percent
- The limit of about 2.5 times for SMT and clock frequency is an estimate, not a measurement

## Re-measurement on 2026-10-24

On a host with 12 CPUs under heavy load from other jobs (load average 37 to 70), user CPU time at 1 and at 12 threads: `rsv/a/2000` 16.3 s and 17.4 s; `dengue/2000` with `--branch-length-mode marginal` 20.1 s and 23.0 s (system time 1.9 s and 5.0 s, wall time 22.1 s and 3.8 s). The ratio of 1.1 to 1.4 is far below the 5.5 to 9.5 above, and the code of the FFT plan cache is unchanged, so no commit explains the change. Spinning threads count for less CPU time on a loaded host, so the issue stays open until a measurement on an idle host with the method above confirms the ratio. The default run of `dengue/2000` now takes about 5 s instead of 120 s, because `--branch-length-mode auto` (9bc74ac83) resolves to input mode on its tree.

## FFT plan cache

The time-grid convolutions call `fn convolve_fft()` in [packages/treetime-ops/src/convolution.rs#L71-L84](../../packages/treetime-ops/src/convolution.rs#L71-L84) (from `fn convolution_function_function_fine()` in [packages/treetime-distribution/src/distribution_ops/convolve.rs#L200](../../packages/treetime-distribution/src/distribution_ops/convolve.rs#L200)):

- The FFT processors of `ndarray-conv` are cached per thread in `static FFT_PROCESSORS_BY_FFT_LEN` [packages/treetime-ops/src/convolution.rs#L10-L15](../../packages/treetime-ops/src/convolution.rs#L10-L15), keyed by transform length, with at most 8 lengths per thread (`FFT_PROCESSORS_PER_THREAD`). Each processor holds its own `realfft` and `rustfft` planners, so each thread plans each transform length that it meets, and a thread that meets more than 8 lengths evicts and plans again
- Transforms longer than 2^14 points (`MAX_CACHED_FFT_LEN`) bypass the cache: `conv_fft()` of `ndarray-conv` builds a new processor, and with it new plans, on every call

Suspected cause: FFT plans built again per thread and per call. The `sincos` share fits the twiddle factors that planning computes (inferred).

> [!IMPORTANT]
> **Investigation required.** The cause of the extra CPU time is not measured. Collect: the number of FFT plans created per run and per thread at 1 and 20 threads, with the transform lengths; CPU profiles of the same case at 1 and 20 threads, compared symbol by symbol; and the parallel efficiency (CPU time and wall time) of each timetree pass. If plan creation is the main cost, share plans between threads (`rustfft` plans are `Arc<dyn Fft>` and can be shared) or widen the cache.
