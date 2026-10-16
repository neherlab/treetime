# Zero branch length boundary is handled in several places

> [!IMPORTANT]
> **Decision required.** The per-edge optimizer can be restructured as "maximize on the closed interval, then compare the interior candidate with the endpoint $t = 0$". This reproduces current Brent outputs only together with a decision on the global search for non-unimodal models, and it reproduces current Newton outputs only after [M-optimize-newton-returns-unverified-point.md](M-optimize-newton-returns-unverified-point.md) is fixed. Options are below.

## Problem

`fn optimize_edge()` [packages/treetime/src/optimize/dispatch.rs#L221-L249](../../packages/treetime/src/optimize/dispatch.rs#L221-L249) computes once per edge whether $t = 0$ is admissible (`ZeroLength::of_edge()` [packages/treetime/src/optimize/zero_boundary.rs#L20-L33](../../packages/treetime/src/optimize/zero_boundary.rs#L20-L33): every site has a positive finite likelihood at zero, and the indels do not require a positive length) and passes the value on. The endpoint $t = 0$ is then handled in these places:

1. **Start adjustment** `fn start_branch_length()` [zero_boundary.rs#L35-L52](../../packages/treetime/src/optimize/zero_boundary.rs#L35-L52): a zero start moves to $\max(k / \text{rate}, \text{one\_mutation})$ when the edge has $k > 0$ indels, or to `one_mutation` when some site has zero likelihood at $t = 0$
2. **Unimodal shortcut** `fn is_zero_branch_optimal()` [zero_boundary.rs#L62-L77](../../packages/treetime/src/optimize/zero_boundary.rs#L62-L77) before the solver: when every partition has the same certified decay rate (`certified_decay_rate()`, JC69, F81, and binary models), zero is admissible, and the edge has no indels, return zero when $\ell'(0) < 0$
3. **Solver floors**: the Brent bracket lower end `min_branch_length.max(1e-12)` [packages/treetime/src/optimize/method_brent.rs#L39](../../packages/treetime/src/optimize/method_brent.rs#L39), so Brent cannot return zero. Commit `fa55c9d6` chose "any positive value below the grid lower bound" and records no further reason: argmin accepts a zero lower end and starts at an interior point, so the floor is not an argmin requirement, and zero needs exclusion only where the objective is undefined at $t = 0$. `min_branch_length_for_indels()` [zero_boundary.rs#L54-L60](../../packages/treetime/src/optimize/zero_boundary.rs#L54-L60) gives $0.01 \cdot \text{one\_mutation}$ for edges that need a positive length. Newton clamps $s \geq 0$ and can return exactly zero
4. **Grid inside Newton** [packages/treetime/src/optimize/method_newton.rs#L28-L30](../../packages/treetime/src/optimize/method_newton.rs#L28-L30): the grid search compares its best positive point with zero
5. **Reconcile** `fn reconcile_zero_boundary()` [zero_boundary.rs#L79-L99](../../packages/treetime/src/optimize/zero_boundary.rs#L79-L99) after the solver: runs the grid search when zero beats a positive candidate, when a zero candidate comes from a model without a certified decay rate, or when a zero candidate is excluded because a site is invalid at zero

## Evidence

Throwaway instrumentation of `fn optimize_edge()` on `rust` at `77dd3a84` (method in [M-optimize-newton-returns-unverified-point.md](M-optimize-newton-returns-unverified-point.md#evidence)). Measured: `treetime optimize` on `flu/h3n2/{20,200}`, `ebola/{20,100}`, `zika/86`, `rsv/a/100`, `--dense true|false`, `--model infer|jc69`, both methods (48 runs, 19,086 per-edge optimizations).

Which path decided the result:

| Path                                |             Brent |           Newton |
| ----------------------------------- | ----------------: | ---------------: |
| Unimodal shortcut returned zero     |               745 |              745 |
| Reconcile changed the solver output | 752 (all to zero) | 2 (both to zero) |
| Solver output used unchanged        |             8,045 |            8,797 |

- **Brent and the reconcile step are coupled.** Every reconcile change in the Brent runs is a solver result within $10^{-4}$ in $s$ of the bracket lower end $s = 10^{-6}$ (at full precision $t$ between $1.5\cdot10^{-10}$ and $4.4\cdot10^{-10}$, so inside the solver tolerance of the endpoint, not on it); with `--model jc69` the same edges take the unimodal shortcut instead (`flu/h3n2/200` dense: 225 such results with `infer`, 225 shortcuts with `jc69`). The positive floor produces the case that the reconcile step then repairs with a 100-point grid search. Whether a zero lower end would remove the case is unmeasured (see the investigation block below)
- **The final result is zero on 1,497 Brent and 1,496 Newton edges.** No edge returned zero while the reference found a better positive point, except 4 Newton edges where the grid inside Newton cannot represent a maximum below $0.1 \cdot \text{one\_mutation}$ ([M-optimize-newton-returns-unverified-point.md](M-optimize-newton-returns-unverified-point.md))
- **On the Dinh and Matsen K80 test surface** (local maximum $t \approx 0.167$, global maximum $t \approx 6.37$), the reconcile step rescues a Newton curvature-break result $t = 1.05$ that is worse than the start: $\ell(0) > \ell(1.05)$, so it runs the grid search, which returns $t \approx 0.162$. An endpoint-only comparison returns zero there. With Brent, the closed-interval design and the current code return the same point for every start in $\{0.001, \dots, 5\}$

> [!IMPORTANT]
> **Investigation required.** The closed-interval design was compared with the current code only on the Dinh and Matsen surface; the dataset comparison stopped before the instrumented build ran. Measure, on the datasets listed above plus `mpox/clade-ii/100`, `dengue/100`, and fixed `k80`/`hky85` models, how many final branch lengths and log-likelihoods change when (a) the Brent bracket starts at $s = 0$ instead of $s = 10^{-6}$ and (b) the reconcile grid search is replaced by the endpoint comparison $\ell(0) \geq \ell(\text{candidate})$. Also count how often the reconcile grid search runs per optimize iteration, to state its cost.

## Options

Both options can change outputs and need approval.

1. **Closed interval plus endpoint comparison**: keep the unimodal shortcut, which is exact; run the solver on $[s_\text{min}, s_\text{max}]$ with $s_\text{min} = 0$ when `ZeroLength` is admissible; return zero when $\ell(0) \geq \ell(\text{candidate})$. Removes the $10^{-12}$ floor and the three reconcile triggers. Leaves non-unimodal models with a local solver only, which loses the grid check that the reconcile step runs today when a candidate is worse than zero: an endpoint comparison cannot replace that search, and positive solutions that the grid finds today can be lost
2. **Option 1 plus an explicit global check for non-unimodal models**: run a grid or multistart search for every non-unimodal edge, independent of the zero endpoint. Separates "find the global maximum" from "handle the boundary", and also covers the case where the global maximum is far from zero (the Dinh and Matsen surface). Costs a grid search per edge for GTR models; its extent is the open question of [N-optimize-reconcile-grid-extent-unverified.md](N-optimize-reconcile-grid-extent-unverified.md). A finite grid or multistart search remains a heuristic, not a global guarantee, so this option needs an approved search domain and accuracy contract

Both options keep the current outputs' differences from v0 open: v0 uses a different bracket and failure policy, see [M-optimize-brent-setup-differs-from-v0.md](M-optimize-brent-setup-differs-from-v0.md#axis-1-search-interval-and-start).

## Cross-references

- [kb/reports/zero-branch-length-optimization.md](../reports/zero-branch-length-optimization.md): description of the current paths
- [N-optimize-reconcile-grid-extent-unverified.md](N-optimize-reconcile-grid-extent-unverified.md), [N-optimize-grid-search-resolution-unverified.md](N-optimize-grid-search-resolution-unverified.md)
- [M-optimize-brent-setup-differs-from-v0.md](M-optimize-brent-setup-differs-from-v0.md)
