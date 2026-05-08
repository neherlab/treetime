# optimize_tree_marginal uses signed convergence check

## v0 location

`TreeAnc.optimize_tree_marginal()` (`#TreeAnc`, `#optimize_tree_marginal`) [packages/legacy/treetime/treetime/treeanc.py#L1297-L1360](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360)

## Erratum

The convergence check at [packages/legacy/treetime/treetime/treeanc.py#L1357](../../packages/legacy/treetime/treetime/treeanc.py#L1357) uses a signed comparison:

```python
deltaLH = LH - oldLH
if deltaLH < LHtol:
    break
```

This accepts two distinct conditions as "convergence":

1. $0 \le \Delta\text{LH} < \text{LHtol}$: small improvement (true convergence)
2. $\Delta\text{LH} < 0$: likelihood decreased (optimization failure)

In soft-EM (<a id="cite-1"></a>[Dempster, Laird, and Rubin 1977](https://doi.org/10.1111/j.2517-6161.1977.tb01600.x) [[1](#ref-1)]), the log-likelihood is guaranteed non-decreasing per iteration. A negative $\Delta\text{LH}$ indicates a violation of this guarantee, which should be diagnosed rather than accepted as convergence. The signed check conflates convergence with failure, reporting both as success (`return ttconf.SUCCESS`).

The correct convergence criterion for an EM-like algorithm is: the likelihood improvement is positive but below a threshold ($0 \le \Delta\text{LH} < \text{LHtol}$). A separate check should handle likelihood decrease (revert to previous state, or report the violation).

## Why v0 is unaffected in practice

v0's `optimize_tree_marginal` operates in a regime where likelihood decreases are vanishingly rare:

- Dense-only mode: v0 has no sparse representation. Dense marginal reconstruction uses full probability matrices (soft-EM), which guarantees monotone likelihood increase before damping.
- Damping breaks monotonicity formally but not practically: exponential damping blends the ML-optimal branch lengths with previous values. This is a relaxation step, not a maximization, so the formal EM guarantee does not hold. In practice, the damped step is a convex combination of the current and previous ML optima, and the likelihood change is almost always positive.
- No in-loop topology changes: v0 calls `prune_short_branches()` after the loop ([packages/legacy/treetime/treetime/treeanc.py#L1434-L1436](../../packages/legacy/treetime/treetime/treeanc.py#L1434-L1436)), not inside it. There are no sudden likelihood drops from edge collapse or sibling merging during iteration.
- No indel contribution: v0 ignores indels in the likelihood. There is no per-iteration indel rate recomputation and no Poisson feedback loop.
- Low iteration count: `max_iter=10` with `damping=0.75` keeps the damping weight at $\ge 0.056$ ($0.75^{10}$) throughout. The optimizer never reaches the fully undamped regime.

Under these conditions, $\Delta\text{LH} < 0$ does not occur in normal operation. The signed check fires exclusively via condition 1 (small positive improvement). The deficient condition 2 (negative delta) is dead code.

## Evidence

- `optimize_tree_marginal` at [packages/legacy/treetime/treetime/treeanc.py#L1357](../../packages/legacy/treetime/treetime/treeanc.py#L1357) uses `deltaLH < LHtol` (signed)
- The adjacent `optimize_tree` (joint mode) at [packages/legacy/treetime/treetime/treeanc.py#L1466](../../packages/legacy/treetime/treetime/treeanc.py#L1466) uses `N_diff < 1` (number of changed nucleotides), a different and more principled criterion that directly measures convergence without confusing it with failure
- The EM convergence theorem (<a id="cite-2"></a>[Wu 1983](https://projecteuclid.org/journals/annals-of-statistics/volume-11/issue-1/On-the-Convergence-Properties-of-the-EM-Algorithm/10.1214/aos/1176346060.full) [[2](#ref-2)], Theorem 2) proves monotone likelihood increase for the standard EM algorithm. The damped variant does not satisfy the theorem's conditions (the M-step does not maximize the Q-function), but the practical monotonicity observation in v0 is consistent with the analysis by <a id="cite-3"></a>[Polyak and Tremba 2020](https://arxiv.org/abs/1703.07810) [[3](#ref-3)] that damped Newton retains monotone progress when the damping factor is close to 1.
- The log message at [packages/legacy/treetime/treetime/treeanc.py#L1358](../../packages/legacy/treetime/treetime/treeanc.py#L1358) prints `deltaLH=%f, stopping iteration` without distinguishing convergence ($\Delta\text{LH} \ge 0$) from failure ($\Delta\text{LH} < 0$)

## v0 impact

- Zero impact under normal conditions (dense mode, damped iterations, no in-loop topology changes)
- Masks any future v0 change that introduces likelihood decreases (e.g. adding sparse mode, in-loop pruning, or indel contributions)

## v1 status

Implemented. `run_optimize_loop()` [packages/treetime/src/commands/optimize/run.rs#L260](../../packages/treetime/src/commands/optimize/run.rs#L260) uses a three-condition check via `ConvergenceReason` enum: converged (small absolute change), oscillating (small 2-step change detecting 2-cycles), or worsened (revert to best-observed branch lengths and stop). See M-optimize-sparse-em-2-cycle (resolved).

## References

- <a id="ref-1"></a>Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. "Maximum Likelihood from Incomplete Data via the EM Algorithm." _Journal of the Royal Statistical Society: Series B_ 39 (1): 1-38 (with discussion). https://doi.org/10.1111/j.2517-6161.1977.tb01600.x [↩](#cite-1)
- <a id="ref-2"></a>Wu, C. F. Jeff. 1983. "On the Convergence Properties of the EM Algorithm." _The Annals of Statistics_ 11 (1): 95-103. https://projecteuclid.org/journals/annals-of-statistics/volume-11/issue-1/On-the-Convergence-Properties-of-the-EM-Algorithm/10.1214/aos/1176346060.full (DOI: https://doi.org/10.1214/aos/1176346060) [↩](#cite-2)
- <a id="ref-3"></a>Polyak, Boris, and Andrey Tremba. 2020. "New Versions of Newton Method: Step-Size Choice, Convergence Domain and Under-Determined Equations." _Optimization Methods and Software_ 35 (6): 1272-1303. https://arxiv.org/abs/1703.07810 (DOI: https://doi.org/10.1080/10556788.2019.1669154) [↩](#cite-3)
