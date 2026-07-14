# Clock SVG output is nondeterministic for an untraced reason

Repeated runs of the same binary produced different root-to-tip SVG bytes. The clock CSV has a confirmed scheduling-dependent row-order defect, but the SVG path has not been traced far enough to determine whether its difference comes from input ordering, floating-point reduction order, metadata, or rendering.

The difference is reproduced at one versus sixteen Rayon threads in [kb/reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-parallel-sparse-leaf-setup.md](../reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-parallel-sparse-leaf-setup.md#output-equivalence) and [kb/reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-zstd.md](../reports/2026-07-13_perf-report-parallel-scaling/mpox-2000-zstd.md#output-equivalence), while duplicate one-thread runs are byte-identical.

## Potential solutions

- O1. Trace both SVGs to their first semantic difference and stabilize the upstream value or ordering responsible for it.
- O2. Normalize the SVG after rendering. This can hide a real numerical or ordering difference and is invalid without proof that only non-semantic metadata varies.

## Recommendation

Use O1. Compare parsed SVG structure and every numeric series before changing code. Create an implementation ticket only after the responsible source value or renderer field is identified.

## Related issues

- [M-clock-parallel-output-order-nondeterministic.md](M-clock-parallel-output-order-nondeterministic.md)
