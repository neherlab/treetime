# Branch length optimization oscillates without damping

The optimize command's outer loop alternates between marginal reconstruction and per-edge branch length optimization without damping. Empirically, this causes oscillation that impedes convergence.

## Problem

The convergence loop in `run_optimize()` ([packages/treetime/src/commands/optimize/run.rs#L121-L141](../../packages/treetime/src/commands/optimize/run.rs#L121-L141)) alternates two steps:

1. Fix branch lengths, recompute marginal ancestral profiles (`update_marginal()`)
2. Fix profiles, optimize all branch lengths via Newton/grid (`run_optimize_mixed()`)

Each step is optimal given the other's output, but the alternation can oscillate: changing all branch lengths shifts marginal profiles, which shifts optimal branch lengths, and so on. Without damping, the full Newton-optimal branch length replaces the previous value each iteration, amplifying oscillations.

## v0 solution: exponential damping

v0 `optimize_tree_marginal()` ([packages/legacy/treetime/treetime/treeanc.py#L1297-L1360](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1360)) applies per-iteration damping:

```python
update_val = new_val * (1 - damping ** (i + 1)) + old_val * damping ** (i + 1)
```

With `damping=0.75` (default), the blend factors are:

| Iteration | New weight | Old weight |
| --------- | ---------- | ---------- |
| 1         | 0.250      | 0.750      |
| 2         | 0.438      | 0.562      |
| 3         | 0.578      | 0.422      |
| 5         | 0.763      | 0.237      |
| 10        | 0.944      | 0.056      |

Early iterations take conservative steps. As the optimization progresses and the landscape stabilizes, damping decays and steps approach the full Newton update. v0 also uses progressive tolerance tightening for Brent: `tol = 1e-8 + 0.01^(i+1)`.

## v1 current state

- No damping at the outer loop level ([packages/treetime/src/commands/optimize/run.rs#L121-L141](../../packages/treetime/src/commands/optimize/run.rs#L121-L141))
- Per-edge Newton step clamps to `[-1.0, current_bl]` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L213](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L213)) which prevents negative branch lengths but does not prevent oscillation
- Grid search fallback (100 points) activates when Newton's Hessian is non-negative, but this does not address alternation oscillation
- `--dp` convergence threshold (default 1e-2) triggers early stop but cannot distinguish converging from oscillating behavior

## Scientific background

Branch length optimization via alternating coordinate maximization (fixing other parameters, optimizing one branch at a time, then re-estimating ancestral states) is standard in phylogenetic ML inference. Clancy, Lyu & Roch (arXiv 2507.22038) prove exponential convergence in favorable regimes with sufficiently close initialization. Their result assumes a strong signal-to-noise regime and does not cover the weak-signal case where oscillation is most likely.

Damped updates are the standard remedy for oscillation in alternating optimization. The damping schedule does not change the fixed point (the full update eventually dominates) but smooths the path toward it.

For the per-edge Newton step specifically, Polyak & Tremba (arXiv 1703.07810) show that adaptively choosing between pure Newton (fast local convergence) and damped Newton (larger convergence domain) gives global convergence guarantees.

## Proposed solutions

### S1: Outer-loop damping (match v0)

Add v0-style exponential damping to the convergence loop. After `run_optimize_mixed()` updates each edge's branch length, blend with the previous value using `damping^(i+1)`. Expose `--damping` CLI parameter (default 0.75).

### S2: Backtracking line search on Newton step

In `run_optimize_mixed()` ([packages/treetime/src/commands/optimize/optimize_unified.rs#L211-L231](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L211-L231)), after computing the Newton step, verify the total log-likelihood actually improved. If not, halve the step and retry. This prevents individual edge overshoot from propagating through the alternation. Can combine with S1.

### S3: Monitor oscillation

Track log-likelihood history. If `lh[i] < lh[i-1]` (regression), reduce step sizes. If the absolute change oscillates sign for >2 consecutive iterations, increase damping automatically.

S1 is the minimal change matching v0 behavior.

## Affected commands

- `optimize` (directly)
- `timetree` (reuses `PartitionOptimizeOps` trait but has its own iteration loop, may need independent damping)
