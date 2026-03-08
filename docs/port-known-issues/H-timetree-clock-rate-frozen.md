# Clock rate frozen across iterations

The clock rate does not update during timetree iterations regardless of
coalescent, marginal, or other settings. The rate is determined by the initial
root-to-tip regression and stays constant. v0 iteratively refines the rate each
iteration.

## Root cause

The clock regression at
[`clock_regression.rs#L74`](../../packages/treetime/src/commands/clock/clock_regression.rs#L74)
reads `edge.branch_length()`, which is the original input branch length from the
Newick file. This value never changes across iterations.

v0 uses `node.clock_length` (time-domain branch length updated by the solver
each iteration). v1 has two separate fields: `edge.branch_length()` (input,
frozen) and `edge.time_length()` (solver output, updated). The regression reads
only the frozen input field.

## Evidence

Verbose log for flu/h3n2/20:

```
Clock rate: 2.842844e-3   (first estimate, before rerooting)
Clock rate: 2.816881e-3   (second estimate, after first ancestral reconstruction)
Clock rate: 2.816881e-3   (all subsequent iterations - never changes)
```

v0 trace shows 3 iterations with varying seq_LH (-3163, -3168, -3161). v1
converges after 1 iteration with n_diff=0.

The rate 0.002817 is produced for all parameter settings on flu/h3n2/20 (default,
coalescent=10, coal-opt, marginal, confidence). v0 produces different rates per
setting (0.002921, 0.002668, 0.002700, 0.002785, 0.002754).

## Impact on clock rate accuracy

Rate divergence from v0 increases with dataset size because larger trees benefit
more from iterative refinement:

| Dataset      | v1 rate  | v0 rate  | Diff   |
| ------------ | -------- | -------- | ------ |
| flu/h3n2/20  | 0.002817 | 0.002921 | -3.6%  |
| flu/h3n2/200 | 0.003354 | 0.003061 | +9.6%  |
| flu/h3n2/500 | 0.003731 | 0.003157 | +18.2% |
| ebola/20     | 0.000626 | 0.000559 | +12.0% |
| zika/20      | 0.001032 | 0.001014 | +1.8%  |

## Impact on root date

The frozen rate compounds into root date discrepancy. flu/h3n2/20 default: v1
root = 1993.47, v0 root = 1997.11 (3.64 years earlier). The earliest tip is
~2000, making the v1 root implausibly early.

| Dataset      | Root v1 | Root v0 | Root diff   |
| ------------ | ------- | ------- | ----------- |
| flu/h3n2/20  | 1993.47 | 1997.11 | -3.64 years |
| flu/h3n2/200 | 1963.44 | 1964.33 | -0.89 years |
| ebola/20     | 2013.01 | 2013.55 | -0.54 years |
| zika/20      | 2012.07 | 2013.00 | -0.93 years |

## Fix direction

Change the clock regression to read `edge.time_length()` (solver output) instead
of `edge.branch_length()` (input). The regression should use the current best
estimate of time-domain branch lengths, not the raw input.

## Related issues

- [Positional likelihood metric differs from v0](M-timetree-positional-likelihood-metric.md)
  interacts with convergence detection
- [Internal node dates missing at scale](M-timetree-internal-dates-missing-scale.md)
  worsened by inaccurate rate at scale
