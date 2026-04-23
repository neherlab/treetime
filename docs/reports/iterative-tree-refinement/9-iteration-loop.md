# Chapter 9: The iteration loop -- putting it all together

[Back to index](_index.md) | Previous: [Chapter 8: Initial estimation](8-initial-estimation.md) | Next: [Chapter 10: Implementation status](10-status-and-recommendations.md)

## The alternating optimization framework

Tree refinement alternates between two interdependent operations:

- **E-step:** Given branch lengths, compute ancestral state distributions via marginal reconstruction ([Chapter 4](4-ancestral-reconstruction.md)). This is Felsenstein's pruning algorithm running backward and forward.
- **M-step:** Given ancestral state distributions, optimize each branch length via Newton-Raphson or Brent's method ([Chapter 5](5-branch-length-optimization.md)).

Undamped EM guarantees monotone likelihood increase per iteration (Dempster, Laird, and Rubin 1977). The sequence converges to a stationary point of the likelihood surface (Wu 1983). Meng and Rubin (1993) proved that updating parameter blocks sequentially (ECM -- Expectation Conditional Maximization) preserves EM's convergence guarantees.

The v1 loop uses heuristic alternating optimization with damping and an absolute delta-LH stop rule. Damping stabilizes convergence empirically but does not carry the formal EM monotonicity guarantee.

## Why it oscillates

Changing all branch lengths in one sweep shifts the ancestral distributions, which shifts the optimal branch lengths in the next sweep. The optimization bounces between two nearby configurations. This is the same phenomenon Young (1954) characterized for the Gauss-Seidel method in linear systems.

## Damping

**Exponential damping** blends new branch lengths with old ones after each sweep:

```
f = max(d^(i+1), DAMPING_FLOOR)
bl_damped = bl_new * (1 - f) + bl_old * f
```

where `d` is the damping factor (default 0.75), `DAMPING_FLOOR = 0.01`, and `i` is the 0-based iteration index. The floor prevents fully undamped late iterations.

| Iteration | New weight | Old weight |
| --------- | ---------- | ---------- |
| 0         | 0.250      | 0.750      |
| 1         | 0.438      | 0.562      |
| 2         | 0.578      | 0.422      |
| 5         | 0.822      | 0.178      |
| 10        | 0.958      | 0.042      |

Early iterations take conservative steps. Later iterations approach the full Newton update. The old-value weight decays exponentially, ensuring convergence to the same fixed point as undamped optimization.

Damping is analogous to **under-relaxation** in SOR (Successive Over-Relaxation, Young 1954). The analogy is informal -- SOR theory applies to linear systems, while phylogenetic optimization is nonlinear. The default `d = 0.75` is empirically chosen, not derived from spectral analysis. The damping factor `d = 0.75` corresponds to relaxation parameter `omega = 0.25`.

## The five loops

### Loop 1: v0 `optimize_tree` (joint mode)

```
Initial: reconstruct ancestral sequences + infer GTR
Loop (max_iter=5):
  1. prune_short_branches()           <-- topology cleanup inside loop
  2. reconstruct_anc()                <-- E-step
  3. if N_diff < 1: break
  4. optimize_branch_lengths_joint()  <-- M-step (no damping)
```

v0 code: [`packages/legacy/treetime/treetime/treeanc.py#L1384-L1473`](../../../packages/legacy/treetime/treetime/treeanc.py#L1384-L1473)

### Loop 2: v0 `optimize_tree_marginal`

```
Loop (max_iter=10, damping=0.75):
  1. for each node: compute damped branch length   <-- M-step with inline damping
  2. infer_ancestral_sequences()                    <-- E-step
  3. if delta_LH < LHtol: break
Then:
  prune_short_branches()                            <-- topology cleanup after loop
```

v0 code: [`packages/legacy/treetime/treetime/treeanc.py#L1297-L1360`](../../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360)

### Loop 3: v0 `_run` (timetree)

```
Loop (max_iter):
  1. add_coalescent_model(Tc)
  2. relaxed_clock()
  3. resolve_polytomies(stochastic_resolve)
  4. if polytomies resolved:
       prepare_tree()
       optimize_tree(max_iter=0, prune_short=False)  <-- re-optimize, no pruning
       make_time_tree()
  5. else:
       infer_ancestral_sequences()
       make_time_tree()
  6. converge when ndiff == 0 AND n_resolved == 0
```

After polytomy resolution, v0 re-optimizes branch lengths but with `prune_short=False` to avoid collapsing newly introduced zero-mutation branches.

v0 code: [`packages/legacy/treetime/treetime/treetime.py#L307-L376`](../../../packages/legacy/treetime/treetime/treetime.py#L307-L376)

### Loop 4: v1 `run_optimize`

```
Initial: initial_guess_mixed()
Loop (max_iter=20):
  1. update_marginal()                            <-- E-step (sparse + dense)
  2. compute total_lh
  3. if delta_lh < dp: break
  4. save_branch_lengths()
  5. run_optimize_mixed()                          <-- M-step (Newton/grid per edge)
  6. find_zero_optimal_internal_edges()            <-- record zero-optimal edges
  7. apply_damping(0.75)                           <-- post-pass damping
  8. prune_and_merge_in_loop(zero_optimal_edges)   <-- topology cleanup
```

Zero-optimal edges are identified after the M-step but before damping (step 6), because damping blends the zero values with old branch lengths and obscures the optimizer's decision. Damping applies normally to all edges (step 7). Then `prune_and_merge_in_loop` overrides the damped values for zero-optimal edges back to zero, collapses them, and merges shared mutations in resulting polytomies (sparse partitions only).

v1 code: [`packages/treetime/src/commands/optimize/run.rs#L275-L378`](../../../packages/treetime/src/commands/optimize/run.rs#L275-L378)

### Loop 5: v1 timetree refinement

```
Loop (convergence-controlled):
  1. apply_relaxed_clock()
  2. capture_ancestral_states()
  3. resolve_polytomies(zero_branch_slope)    <-- greedy temporal
  4. prepare_tree_after_topology_change()
  5. run_timetree() + update_marginal()
  6. count_sequence_changes()
  7. re-estimate clock model
```

v1 code: [`packages/treetime/src/commands/timetree/refinement.rs#L21-L106`](../../../packages/treetime/src/commands/timetree/refinement.rs#L21-L106)

## Cross-loop comparison

| Aspect                       | v0 joint    | v0 marginal     | v0 timetree                 | v1 optimize        | v1 timetree         |
| ---------------------------- | ----------- | --------------- | --------------------------- | ------------------ | ------------------- |
| Prune zero-length            | Inside loop | After loop      | After resolution (no prune) | Inside loop        | Not implemented     |
| Polytomy resolution          | N/A         | N/A             | Greedy or stochastic        | Merge after prune  | Greedy only         |
| Damping                      | None        | Inline per-edge | None                        | Post-pass blending | N/A                 |
| Convergence                  | N_diff < 1  | delta_LH < tol  | ndiff==0 AND n_resolved==0  | delta_LH < dp      | n_diff + n_resolved |
| Branch re-opt after topology | N/A         | N/A             | Yes                         | N/A                | No                  |

## Implementation notes

The topology cleanup in Loop 4 identifies zero-optimal edges before damping (step 6) rather than after. This avoids using damped branch lengths as the pruning criterion: damping blends zero with the old value, so a damped zero-optimal edge has a small positive branch length that does not reflect the optimizer's actual decision. The recorded zero-optimal edge keys are used after damping to override the damped values back to zero and collapse the edges.

## Convergence theory

The undamped EM guarantee (Dempster, Laird, and Rubin 1977; Wu 1983) ensures monotonic likelihood increase per iteration. Tseng (2001) extends this to block coordinate descent. Wright (2015) surveys convergence rates. The v1 loop uses damping and an absolute delta-LH stop rule, which stabilizes convergence empirically but does not carry the formal EM monotonicity guarantee.

The **resolution threshold** (0.05 for polytomy resolution) and the **pruning criterion** (branch length < 0.1/L AND zero-distance probability > 0.1) are heuristic safeguards against the star tree paradox (Steel and Penny 2000) -- they prevent the optimizer from endlessly creating and collapsing structure at the resolution boundary.

## References

- Dempster, A. P., N. M. Laird, and D. B. Rubin. 1977. "Maximum Likelihood from Incomplete Data Via the EM Algorithm." _JRSS:B_ 39(1):1-38. https://doi.org/10.1111/j.2517-6161.1977.tb01600.x
- Wu, C. F. J. 1983. "On the Convergence Properties of the EM Algorithm." _Ann. Stat._ 11(1):95-103. https://projecteuclid.org/journals/annals-of-statistics/volume-11/issue-1/On-the-Convergence-Properties-of-the-EM-Algorithm/10.1214/aos/1176346060.full (DOI: https://doi.org/10.1214/aos/1176346060)
- Meng, X.-L., and D. B. Rubin. 1993. "Maximum Likelihood Estimation via the ECM Algorithm." _Biometrika_ 80(2):267-278. https://doi.org/10.1093/biomet/80.2.267
- Young, D. M. 1954. "Iterative Methods for Solving Partial Difference Equations." _Trans. Amer. Math. Soc._ 76(1):92-111. https://doi.org/10.1090/S0002-9947-1954-0059635-7
- Tseng, P. 2001. "Convergence of a Block Coordinate Descent Method." _J. Optim. Theory Appl._ 109(3):475-494. https://doi.org/10.1023/A:1017501703105
- Wright, S. J. 2015. "Coordinate Descent Algorithms." _Math. Program._ 151(1):3-34. https://doi.org/10.1007/s10107-015-0892-3
- Steel, M., and D. Penny. 2000. "Parsimony, Likelihood, and the Role of Models." _Mol. Biol. Evol._ 17(6):839-850. https://doi.org/10.1093/oxfordjournals.molbev.a026364
- Sagulenko, P., V. Puller, and R. A. Neher. 2018. "TreeTime." _Virus Evol._ 4(1):vex042. https://doi.org/10.1093/ve/vex042
