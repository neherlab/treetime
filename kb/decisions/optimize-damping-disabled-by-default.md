# The optimize loop is undamped by default

The outer loop of `treetime optimize` alternates marginal reconstruction and per-edge branch length optimization. Damping blends each update with the previous branch length, $b = b_{new}(1-d) + b_{old}\,d$ with $d = \max(\mathrm{damping}^{i+1}, 0.01)$ in iteration $i$. v1 applies no damping by default: `--damping` defaults to 0, which skips the blend and its floor. v0 damps with 0.75.

**Type**: Scientific default that differs from v0. Approved by the maintainer team.

**v0 location**: `optimize_tree_marginal()` takes `damping=0.75` at [`packages/legacy/treetime/treetime/treeanc.py#L1297-L1298`](../../packages/legacy/treetime/treetime/treeanc.py#L1297-L1298) and blends the root arc at [`#L1334`](../../packages/legacy/treetime/treetime/treeanc.py#L1334) and every other branch at [`#L1342`](../../packages/legacy/treetime/treetime/treeanc.py#L1342).

**v1 location**: the default of `damping` in `struct TreetimeOptimizeArgsRaw` in [`packages/app-commands/src/commands/optimize/args.rs`](../../packages/app-commands/src/commands/optimize/args.rs); `fn apply_damping()` and `DAMPING_FLOOR` in [`packages/treetime/src/optimize/iteration.rs`](../../packages/treetime/src/optimize/iteration.rs).

## Rationale

Damping suppressed the 2-cycle of edges in series: two edges joined by a node without data enter the likelihood only through their sum, and optimizing them independently maps the sum $s$ to $2s^* - s$. v1 optimizes the sum of each such chain once per iteration ([optimize-series-edges-joint-sum.md](optimize-series-edges-joint-sum.md)), so the undamped loop reaches the optimal sum for fixed messages in one step. Damping then only slows convergence. `--damping` stays available for inputs where the undamped loop does not settle.

With damping 0 the floor of 0.01 does not apply either. The loop keeps its protection against a non-monotone likelihood, which the sparse backend can produce: it stops on a converged, oscillating, or worsened likelihood and returns the branch lengths of the iteration with the highest likelihood.

## Evidence

- Simulated data from `validation/` (24 sets, `brent`): damping 0 reaches the same likelihood and the same median branch-length accuracy as 0.75, in 3.3 instead of 8.5 iterations and about half the run time
- Example datasets in `data/` (default settings, sparse and dense on the 20-sample sets, sparse on `lassa/L/50`, `ebola/100`, `flu/h3n2/200`, `mpox/clade-ii/100`, `rsv/a/100`): 78 likelihood evaluations at damping 0 against 115 at 0.75. The highest likelihood of a run is higher at damping 0 in 18 of 21 runs. No run alternates between two likelihoods; every run stops on convergence or after a single decline
- The three runs with a lower likelihood at damping 0 are sparse: `lassa/L/20` by 3.36 (also with 50 iterations and `--tolerance=1e-4`, where both settings stop on a worsened likelihood), `mpox/clade-ii/100` by 0.03, `rsv/a/20` by 0.02. The dense runs of `lassa/L/20` agree within 0.002. `ebola/100` sparse differs by 2.41 in the other direction, so the sparse differences are path dependence of the sparse approximation rather than a systematic loss of the undamped loop
