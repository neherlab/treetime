# Optimize per-branch lengths diverge from v0 fixture beyond 10% relative tolerance

## Problem

The optimize golden-master fixture `packages/treetime/src/optimize/__tests__/__fixtures__/gm_optimize_outputs.json` records a per-branch `final_branch_lengths` map captured from v0 on the same tree, alignment, GTR model, and damping schedule. The summed branch length matches v0 within 5% (asserted by `test_gm_optimize`), but a per-branch comparison after skipping branches v0 collapsed to near-zero (`expected_bl < 1e-5`) exceeds 10% relative tolerance on individual non-zero branches. Example on `flu_h3n2_20_jc69_damped`, the only case of the per-branch test: v1 reports `0.00494` on a branch where v0 reports `0.00071` (a factor of 7 disagreement on a branch larger than the `prune_short_branches` threshold).

## Evidence

`test_gm_optimize_per_branch` in [packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L50-L91](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L50-L91) is `#[ignore]`'d with this issue path in the ignore reason. It runs Brent with `IndelPolicy::Include`. `./dev/docker/run just test-rs --run-ignored only test_gm_optimize_per_branch` reproduces the divergence. The test stops at the first divergent branch, so the full list of divergent branches is not known. CI checks only the summed length within 5% (`test_gm_optimize`, [packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L23-L48](../../packages/treetime/src/optimize/__tests__/test_gm_optimize.rs#L23-L48)).

## Likely contributors

1. Topology cleanup criterion. v0's `prune_short_branches` collapses any internal branch with `bl < 0.1 * one_mutation && prob_t(parent, child, 0) > 0.1`. v1 applies the same rule after the loop (`fn collapse_short_internal_edges()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs)): it collapses internal edges shorter than 0.1/L whose marginal estimate carries no substitution, and every internal edge of length zero. Without this rule, branches that the optimizer leaves at e.g. `1e-12 - 1e-3` subs/site stay positive and absorb mass that v0 redistributes to siblings.

> [!IMPORTANT]
> **Investigation required.** The per-branch divergence has not been measured since the post-loop collapse and the sparse missing-data fix. Run `test_gm_optimize_per_branch` again, list every divergent branch, and check which of the contributors below remain.
2. Indel term. The per-branch test runs with `IndelPolicy::Include`, so v1 adds a Poisson indel term to each per-edge objective. `IndelPolicy::edge_indels()` returns `EdgeIndels` with one indel rate per alignment, estimated once at the start of `fn run_optimize_loop()` ([packages/treetime/src/optimize/run_loop.rs#L47](../../packages/treetime/src/optimize/run_loop.rs#L47)); the loop recounts the indel events per iteration but keeps the rates. v0 has no indel term, so indel-bearing edges get a different objective.
3. Damping schedule application order. v0 and v1 both apply `bl = bl_new * (1 - damping^(i+1)) + bl_old * damping^(i+1)`, but v1 restores zeros for internal zero-optimal edges after damping (see `prune_and_merge_in_loop`), while v0 prunes at the end of each pass with the threshold above. Different cleanup ordering can leave different mass distributions per branch even when the sum agrees.

## Impact

`test_gm_optimize` summed-total check still passes (within 5%), so v1's overall scale is correct. The per-branch disagreement matters for downstream consumers that read individual branch lengths (e.g., comparison runs against v0, ancestral-state visualizations that color edges by length, or downstream timetree inputs that consume the optimized tree).

## Approach

Measure the divergence again with the post-loop collapse in place. If branches still differ, trace them to one of the remaining contributors (indel term, cleanup order). Re-deriving the GM fixture from v1 itself would turn the v0 oracle into a snapshot and needs an approved entry in `decisions/`.

## Cross-references

- [kb/decisions/command-prune-standalone.md](../decisions/command-prune-standalone.md) (v1 prune as a standalone command, not in the optimize loop)
- [`packages/treetime/src/optimize/topology/collapse.rs`](../../packages/treetime/src/optimize/topology/collapse.rs) (shared `collapse_edge()`)
- `fn find_zero_optimal_internal_edges()` and `fn prune_and_merge_in_loop()` [packages/treetime/src/optimize/run_loop.rs#L312-L385](../../packages/treetime/src/optimize/run_loop.rs#L312-L385)
- `def prune_short_branches()` [packages/legacy/treetime/treetime/treeanc.py#L1475-L1495](../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1495)
