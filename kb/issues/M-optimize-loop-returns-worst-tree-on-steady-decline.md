# The optimize loop returns the last tree when the likelihood declines from the start

`fn run_optimize_loop` in [packages/treetime/src/optimize/run_loop.rs](../../packages/treetime/src/optimize/run_loop.rs) reverts to the best branch lengths only when a drop follows a new best (`i >= 2 && lh < lh_prev && lh_prev >= best_lh`). When the likelihood is best at the first iteration and declines steadily afterwards, the condition never holds: the loop runs to `max_iter` and returns the last, lowest-likelihood branch lengths.

## Evidence

Before the branch objective used the correct transition orientation, `optimize` on simulated data (validation/, exp rep07, L=1000, mu=0.01, started from better branch lengths) lowered the log likelihood at every iteration from -11542.8 to -11565.5 and returned the final tree. The objective is fixed, but simultaneous per-branch updates still need not increase the likelihood, so the loop should return the best branch lengths it saw whenever it ends.
