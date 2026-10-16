# Zero branch length boundary is handled in several places

> [!IMPORTANT]
> **Decision required.** The per-edge optimizer can be restructured as "maximize on the closed interval, then compare the interior candidate with the endpoint $t = 0$". This reproduces current Brent outputs only together with a decision on the global search for non-unimodal models, and it reproduces current Newton outputs only after [M-optimize-newton-returns-unverified-point.md](M-optimize-newton-returns-unverified-point.md) is fixed. Options are below.

## Problem

`fn optimize_edge()` [packages/treetime/src/optimize/dispatch.rs#L272-L311](../../packages/treetime/src/optimize/dispatch.rs#L272-L311) handles the endpoint $t = 0$ of the branch length in these places:

1. **Start bumps** [dispatch.rs#L280-L290](../../packages/treetime/src/optimize/dispatch.rs#L280-L290): a zero start moves to $\max(k / \text{rate}, \text{one\_mutation})$ when the edge has $k > 0$ indels, or to `one_mutation` when some site has zero likelihood at $t = 0$
2. **Unimodal shortcut** `fn is_zero_branch_optimal()` [packages/treetime/src/optimize/zero_boundary.rs#L21-L43](../../packages/treetime/src/optimize/zero_boundary.rs#L21-L43) before the solver: for unimodal models (JC69, F81, binary) without indels, return zero when $\ell'(0) < 0$
3. **Solver floors**: the Brent bracket lower end `min_branch_length.max(1e-12)` [packages/treetime/src/optimize/method_brent.rs#L39](../../packages/treetime/src/optimize/method_brent.rs#L39), so Brent cannot return zero; `min_branch_length_for_indels()` [zero_boundary.rs#L13-L19](../../packages/treetime/src/optimize/zero_boundary.rs#L13-L19) gives $0.01 \cdot \text{one\_mutation}$ for edges that need a positive length. Newton clamps $s \geq 0$ and can return exactly zero
4. **Grid inside Newton** [packages/treetime/src/optimize/method_newton.rs#L27-L29](../../packages/treetime/src/optimize/method_newton.rs#L27-L29): the grid search compares its best positive point with zero
5. **Reconcile** `fn reconcile_zero_boundary()` [zero_boundary.rs#L45-L68](../../packages/treetime/src/optimize/zero_boundary.rs#L45-L68) after the solver: runs the grid search when a positive candidate is worse than zero, when a zero candidate comes from a non-unimodal model, or when a zero candidate is invalid at zero

### Repeated scans

`OptimizationContribution::all_sites_valid_at_zero()` [packages/treetime/src/partition/optimize/contribution.rs#L51-L56](../../packages/treetime/src/partition/optimize/contribution.rs#L51-L56) scans every site of every partition (as costly as one log-likelihood evaluation without the exponentials). Per edge it runs in the start bump (zero start only), in `is_zero_branch_optimal` (unimodal models), in `reconcile_zero_boundary` (always), in `is_zero_better_than_grid_best` [zero_boundary.rs#L115](../../packages/treetime/src/optimize/zero_boundary.rs#L115) (positive candidate), and again in `is_zero_better_than_grid_best` at the end of every grid search. The result depends only on the edge's contributions, so one scan per edge is enough.

### Misleading names

- `is_zero_better_than_grid_best(contributions, indels, best_positive)` receives the solver candidate in `reconcile_zero_boundary`, not a grid point
- `branch_length_for_extent` of `reconcile_zero_boundary` is the bumped start, not the input branch length

## Evidence

Throwaway instrumentation of `fn optimize_edge()` on `rust` at `77dd3a84` (method in [M-optimize-newton-returns-unverified-point.md](M-optimize-newton-returns-unverified-point.md#evidence)). Measured: `treetime optimize` on `flu/h3n2/{20,200}`, `ebola/{20,100}`, `zika/86`, `rsv/a/100`, `--dense true|false`, `--model infer|jc69`, both methods (48 runs, 19,086 per-edge optimizations).

Which path decided the result:

| Path                                |             Brent |           Newton |
| ----------------------------------- | ----------------: | ---------------: |
| Unimodal shortcut returned zero     |               745 |              745 |
| Reconcile changed the solver output | 752 (all to zero) | 2 (both to zero) |
| Solver output used unchanged        |             8,045 |            8,797 |

- **Brent and the reconcile step are coupled.** Every reconcile change in the Brent runs is a solver result at the bracket lower end $s = 10^{-6}$; with `--model jc69` the same edges take the unimodal shortcut instead (`flu/h3n2/200` dense: 225 bracket-end results with `infer`, 225 shortcuts with `jc69`). The $10^{-12}$ floor exists only because the argmin bracket is open at zero; it creates the case that the reconcile step then repairs with a 100-point grid search
- **The final result is zero on 1,497 Brent and 1,496 Newton edges.** No edge returned zero while the reference found a better positive point, except 4 Newton edges where the grid inside Newton cannot represent a maximum below $0.1 \cdot \text{one\_mutation}$ ([M-optimize-newton-returns-unverified-point.md](M-optimize-newton-returns-unverified-point.md))
- **On the Dinh and Matsen K80 test surface** (local maximum $t \approx 0.167$, global maximum $t \approx 6.37$), the reconcile step rescues a Newton curvature-break result $t = 1.05$ that is worse than the start: $\ell(0) > \ell(1.05)$, so it runs the grid search, which returns $t \approx 0.162$. An endpoint-only comparison returns zero there. With Brent, the closed-interval design and the current code return the same point for every start in $\{0.001, \dots, 5\}$

> [!IMPORTANT]
> **Investigation required.** The closed-interval design was compared with the current code only on the Dinh and Matsen surface; the dataset comparison stopped before the instrumented build ran. Measure, on the datasets listed above plus `mpox/clade-ii/100`, `dengue/100`, and fixed `k80`/`hky85` models, how many final branch lengths and log-likelihoods change when (a) the Brent bracket starts at $s = 0$ instead of $s = 10^{-6}$ and (b) the reconcile grid search is replaced by the endpoint comparison $\ell(0) \geq \ell(\text{candidate})$. Also count how often the reconcile grid search runs per optimize iteration, to state its cost.

## Options

1. **Closed interval plus endpoint comparison**: compute "zero is admissible" (no required positive length, all sites valid at zero) once per edge; keep the unimodal shortcut, which is exact; run the solver on $[s_\text{min}, s_\text{max}]$ with $s_\text{min} = 0$ when zero is admissible; return zero when $\ell(0) \geq \ell(\text{candidate})$. Removes the $10^{-12}$ floor, the three reconcile triggers, and the repeated scans. Leaves non-unimodal models with a local solver only, which loses the grid check that the reconcile step runs today when a candidate is worse than zero
2. **Option 1 plus an explicit global check for non-unimodal models**: run a grid or multistart search for every non-unimodal edge, independent of the zero endpoint. Separates "find the global maximum" from "handle the boundary", and also covers the case where the global maximum is far from zero (the Dinh and Matsen surface). Costs a grid search per edge for GTR models; its extent is the open question of [N-optimize-reconcile-grid-extent-unverified.md](N-optimize-reconcile-grid-extent-unverified.md)
3. **Keep the structure, remove the repeated scans and rename**: compute `all_sites_valid_at_zero` once and pass it; rename the `best_positive` and `branch_length_for_extent` parameters. No output change

## Cross-references

- [kb/reports/zero-branch-length-optimization.md](../reports/zero-branch-length-optimization.md): description of the current paths (its line links are stale)
- [N-optimize-reconcile-grid-extent-unverified.md](N-optimize-reconcile-grid-extent-unverified.md), [N-optimize-grid-search-resolution-unverified.md](N-optimize-grid-search-resolution-unverified.md)
- [N-optimize-brent-ignores-termination-status.md](N-optimize-brent-ignores-termination-status.md)
