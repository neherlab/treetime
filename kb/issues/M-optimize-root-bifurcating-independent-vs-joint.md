# Bifurcating root optimization: independent per-edge vs joint arc optimization

## Summary

`fn run_optimize_mixed_inner()` optimizes root-incident edges independently in the per-edge loop, then redistributes the total by the pre-optimization ratio. v0 skips root edges from the per-edge loop and jointly optimizes the total root arc length using profiles from both subtrees. These approaches produce different optimized totals.

Two additional issues in the redistribution logic:

1. Post-loop redistribution unconditionally overrides edges the optimizer set to zero via `is_zero_branch_optimal()`.
2. The independent optimization is a form of coordinate descent: it converges to the same stationary point as joint optimization across iterations, but the per-iteration total differs from v0's joint optimum.

## Scientific Background

For time-reversible substitution models (JC69, HKY85, GTR), the likelihood depends only on the total distance between the two subtrees adjacent to the root, not on how that distance is partitioned between the two edges (Felsenstein 1981, "Pulley Principle"). The transition probability satisfies $P(t_1) P(t_2) = P(t_1 + t_2)$ for reversible $Q$, making only the sum $t_1 + t_2$ identifiable from sequence data.

Standard phylogenetic tools (RAxML, IQ-TREE, PhyML) operate on unrooted trees to avoid this ambiguity. When operating on rooted trees, the principled approach is to optimize the root arc jointly as a single parameter.

## v0 Behavior

In `fn optimize_tree_marginal()` [[src](https://github.com/neherlab/treetime/blob/1a5c93d4f7929c3cc9302859015a32562dae02fd/packages/legacy/treetime/treetime/treeanc.py#L1317-L1337)], when the root has exactly two children:

1. Captures `total_bl = n1.branch_length + n2.branch_length` and `bl_ratio = n1.branch_length / total_bl`
2. Constructs profiles from both subtrees: `prof_c = n1.marginal_subtree_LH`, `prof_p = normalize_profile(n2.marginal_subtree_LH * root.marginal_outgroup_LH)`
3. Calls `optimal_t_compressed((prof_p, prof_c), ...)` once for the total arc length
4. Applies damping to the total: `update_val = new_bl * (1 - d^(i+1)) + total_bl * d^(i+1)`
5. Distributes by ratio: `n1.bl = update_val * bl_ratio`, `n2.bl = update_val * (1 - bl_ratio)`

v0 quirk: the loop visits each clade, and both root children satisfy the condition `n.up.up is None and len(n.up.clades) == 2`. The joint optimization block executes twice per iteration (once per child), with updated branch lengths on the second pass but unchanged profiles. This effectively applies double damping to the root arc.

## v1 Behavior

In `fn run_optimize_mixed_inner()` [[src](https://github.com/neherlab/treetime/blob/27ad41454083ddbc09e6ef654063bb3e2b588137/packages/treetime/src/optimize/dispatch.rs#L75-L93)]:

1. Captures ratio from current branch lengths before the loop
2. Optimizes ALL edges (including both root edges) independently in the per-edge loop
3. After the loop, reads the independently optimized branch lengths and redistributes the total by the original ratio [[src](https://github.com/neherlab/treetime/blob/27ad41454083ddbc09e6ef654063bb3e2b588137/packages/treetime/src/optimize/dispatch.rs#L232-L238)]

Each root edge's optimization uses the marginal profiles at both endpoints. The root endpoint's profile incorporates information from the other subtree through the message-passing algorithm. So v1 does not use "partial information" -- it uses coordinate descent (optimize $t_0$ with $t_1$ fixed via profiles, then $t_1$ with $t_0$ fixed) rather than joint optimization (optimize $t_{total}$ directly).

## Issue 1: Zero-branch override

If `is_zero_branch_optimal()` determines one root edge should be zero [[src](https://github.com/neherlab/treetime/blob/27ad41454083ddbc09e6ef654063bb3e2b588137/packages/treetime/src/optimize/dispatch.rs#L135-L138)], the post-loop redistribution unconditionally assigns `total * ratio` (non-zero) to that edge, overriding the optimizer's decision. v0 does not have this problem because root edges are excluded from the per-edge loop.

## Issue 2: Coordinate descent vs joint optimization

The likelihood $L = \sum_r [\pi(r) \cdot P(\text{data}_1 | r, t_1) \cdot P(\text{data}_2 | r, t_2)]$ couples $t_1$ and $t_2$ through the root state summation. Optimizing each independently (coordinate descent) converges to the same stationary point as joint optimization across the outer iteration loop, but the per-iteration total $t_1^* + t_2^*$ may differ from v0's joint optimum $t_{total}^*$. The magnitude depends on tree asymmetry and branch lengths.

## Damping Interaction

The ratio IS preserved through damping. Proof: after redistribution, root edges have ratio `r`. The saved pre-optimization values also have ratio `r` (by construction for iteration 1; by induction for subsequent iterations, since the previous iteration's damping preserves the ratio). Per-edge damping computes `damped = w * redistributed + (1-w) * old`, and since both terms have the same ratio, the result does too.

The damped TOTAL differs from v0's damped total because the independently optimized total differs from v0's jointly optimized total. This is a consequence of Issue 2 above, not an independent problem. Fixing Issue 2 (joint optimization) would make the damped totals match.

Formal: let `ratio = old_0 / (old_0 + old_1)`. After redistribution: `edge_0 = total_opt * ratio`. After damping: `damped_0 = w * total_opt * ratio + (1-w) * old_0 = w * total_opt * ratio + (1-w) * ratio * old_total = ratio * (w * total_opt + (1-w) * old_total)`. Similarly for edge_1 with `(1-ratio)`. The damped ratio equals `ratio`. QED.

## Options

1. Match v0: skip both root-incident edges from the per-edge loop, jointly optimize the total root arc length using contributions from both subtrees, distribute by ratio. Resolves both Issue 1 and Issue 2.
2. Fix Issue 1 only: skip redistribution when either root edge was set to zero. Accept the coordinate descent approach as an intentional change (converges to same point eventually, may take more iterations). Document as a decision.
3. Current approach unchanged: accept both issues as negligible in practice (the difference may be small for well-behaved trees). Quantify with golden master comparison.

## Per-Iteration Divergence

The independent optimization produces a per-iteration root total that is a reflection of the old total about the joint optimum $t^*$:

$$t_0^* + t_1^* = 2t^* - t_{\text{old}}$$

After damping with weight $d$, the per-iteration divergence from v0 is:

$$v1 - v0 = (1 - d) \cdot (t^* - t_{\text{old}})$$

With the default $d = 0.75$, the excess is $0.25 (t^* - t_{\text{old}})$ per step, contracting as the total approaches $t^*$.

At $d = 0$ (damping disabled, a supported CLI option), the map $s_{n+1} = 2t^* - s_n$ is a pure period-2 reflection that never converges. The loop's `Oscillating` stop criterion would halt at an arbitrary phase of the cycle. Joint optimization (Option 1) makes $t^*$ a true fixed point independent of damping.

## Quantification Needed

Empirical comparison of root branch lengths between v0 and v1 on existing datasets (`flu/h3n2/20`, `ebola`, `zika`) would determine whether the coordinate descent approach produces materially different results or converges quickly enough that the per-iteration difference is negligible.
