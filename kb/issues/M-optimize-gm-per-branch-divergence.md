# Optimize per-branch lengths diverge from v0 fixture beyond 10% relative tolerance

## Problem

The optimize golden-master fixture `packages/treetime/src/optimize/__tests__/__fixtures__/gm_optimize_outputs.json` records a per-branch `final_branch_lengths` map captured from v0 on the same tree, alignment, GTR model, and damping schedule. The summed branch length matches v0 within 5% (asserted by `test_gm_optimize`), but a per-branch comparison after skipping branches v0 collapsed to near-zero (`expected_bl < 1e-5`) exceeds 10% relative tolerance on individual non-zero branches. Example on `flu_h3n2_20_jc69_damped`: v1 reports `0.00488` on a branch where v0 reports `0.00071` (a factor of 7 disagreement on a branch larger than the `prune_short_branches` threshold).

## Evidence

`packages/treetime/src/optimize/__tests__/test_gm_optimize.rs` `test_gm_optimize_per_branch` is `#[ignore]`'d with this issue path in the ignore reason. Removing the ignore reproduces the divergence on every non-skipped fixture case.

## Likely contributors

1. Topology cleanup criterion. v0's `prune_short_branches` collapses any
   internal branch with `bl < 0.1 * one_mutation && prob_t(parent, child, 0) > 0.1`. v1 collapses only branches the optimizer drives to exact zero. Branches the optimizer leaves at e.g. `1e-12 - 1e-3` subs/site that v0 collapses stay positive in v1 and absorb mass that v0 redistributes to siblings.
2. Indel rate handling. v1 estimates `indel_rate` once at the start of
   `run_optimize_mixed` and treats it as a constant; v0 may recompute per iteration. A different rate biases the per-edge Newton/Brent objective.
3. Damping schedule application order. v0 and v1 both apply
   `bl = bl_new * (1 - damping^(i+1)) + bl_old * damping^(i+1)`, but v1 restores zeros for internal zero-optimal edges after damping (see `prune_and_merge_in_loop`), while v0 prunes at the end of each pass with the threshold above. Different cleanup ordering can leave different mass distributions per branch even when the sum agrees.

## Impact

`test_gm_optimize` summed-total check still passes (within 5%), so v1's overall scale is correct. The per-branch disagreement matters for downstream consumers that read individual branch lengths (e.g., comparison runs against v0, ancestral-state visualizations that color edges by length, or downstream timetree inputs that consume the optimized tree).

## Approach

Three orthogonal options:

1. Implement v0's `prune_short_branches` criterion (`bl < 0.1 * one_mutation`
   AND `prob_t(parent, child, 0) > 0.1`) inside the v1 optimize loop and recheck the per-branch divergence.
2. Re-derive the GM fixture from v1 itself (snapshot rather than oracle)
   and accept the per-branch divergence as an intentional v0/v1 deviation, moving the entry to `decisions/`.
3. Tighten the v1 internal-edge collapse criterion to match v0's threshold
   (without copying the probability gate).

Option 1 is the most aligned with the porting goal. Option 2 is the fastest path to a green CI but loses the v0 oracle. Option 3 is a compromise but introduces a new threshold parameter.

## Cross-references

- `../decisions/command-prune-standalone.md` (v1 prune as
  a standalone command, not in the optimize loop)
- [`packages/treetime/src/optimize/topology/collapse.rs`](../../packages/treetime/src/optimize/topology/collapse.rs) (shared `collapse_edge()`)
- `packages/treetime/src/commands/optimize/run.rs` `find_zero_optimal_internal_edges`
  and `prune_and_merge_in_loop`
- `packages/legacy/treetime/treetime/treeanc.py` `prune_short_branches` (line 1475)
- Sparse EM 2-cycle non-convergence (fixed): v1 defaults now match v0 (`max_iter=10`, `dp=0.1`), which may have contributed to per-branch divergence
