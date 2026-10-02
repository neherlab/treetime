# Mugration golden master parity with v0

v1 implements iterative GTR inference for mugration, matching v0's `reconstruct_discrete_traits()` ([packages/legacy/treetime/treetime/wrappers.py#L653-L811](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811)) algorithm structure. With the default (v0) inference policy, v1 reproduces v0 trait assignments for zika, zika-with-weights, and lassa. The remaining datasets (dengue, tb, rsv, mpox) still diverge at a few ambiguous internal nodes, and all confidence profiles diverge from v0 by ~1e-3.

## Resolved contributors

- GTR rate-optimizer regression: the rate optimizer briefly used argmin `BrentOpt` (golden-section seeded), which converged to a different rate than v0's interior-seeded scipy `brent` and flipped assignments at ambiguous nodes. Fixed by `BrentBracketed` at [packages/treetime/src/gtr/brent_bracketed.rs](../../packages/treetime/src/gtr/brent_bracketed.rs).
- D1 (initial-pi pseudo-count) and D2 (uninformative-root filtering) are now opt-in flags defaulting to v0 (`--smooth-initial-pi`, `--filter-uninformative-root`). They no longer perturb the default reconstruction.

## Remaining divergence: residual ~1e-3 marginal-confidence difference

The unresolved cause is a ~1e-3 difference in v1's marginal confidence profiles relative to v0, present even for unweighted, informative-root datasets where D1 and D2 are no-ops. Example (zika_20_country root `NODE_0000000`): v0 `0.4812`, v1 `0.4800`. This small difference is below the old ~1.34e-2 figure but still tips the argmax at near-tied nodes, which is why dengue/tb/rsv/mpox assignments differ. The origin (marginal/GTR numerics) is not yet localized.

Per project rules, numerical error > 1e-6 against the v0 oracle is a defect. This remains open.

## Candidate causes (not yet tested against the golden masters)

- **Stale forward messages**: v0 runs a full marginal pass under the first fitted GTR before its rate-and-GTR iterations (`infer_ancestral_sequences` calls `_ml_anc` after `infer_gtr`, [packages/legacy/treetime/treetime/treeanc.py#L564-L566](../../packages/legacy/treetime/treetime/treeanc.py#L564-L566)). v1 `refine_gtr_model_and_rate` ([packages/treetime/src/gtr/refinement.rs](../../packages/treetime/src/gtr/refinement.rs)) refreshes the backward messages through the rate optimization but keeps the forward messages from the initial GTR for every iteration. The claim in [kb/proposals/mugration-full-reconstruction-per-iteration.md](../proposals/mugration-full-reconstruction-per-iteration.md) that v1 matches v0 here does not hold
- **Rate-optimizer evaluation point**: v0's `optimize_gtr_rate` sets the rate to the optimum without a new backward pass, so its subtree messages come from the last Brent evaluation; v1 recomputes them at the optimum

## Leaf evidence

Discrete leaves keep their observed trait as evidence in every pass ([packages/treetime/src/partition/marginal/discrete/partition.rs](../../packages/treetime/src/partition/marginal/discrete/partition.rs)), as v0 does, so a leaf posterior never feeds back as evidence. This is not the cause of the divergence above:

- The golden master datasets have no leaf with a missing trait. Their outputs do not depend on this rule, and zika_20_country root `NODE_0000000` is `0.4800` (v0 `0.4812`)
- Inputs with missing traits diverge more. On `data/zika/20` with the country of three leaves set to `?`, v1 fits the rate 7.7146 and v0 fits 6.4068. The cause is not localized

## Affected golden master tests

- `test_gm_mugration_outputs`: zika, zika_weights, lassa pass; dengue, tb, rsv, mpox ignored (`test_gm_mugration_outputs_v1_divergence`).
- `test_gm_mugration_confidence_zika` (1e-6) and `test_gm_mugration_confidence_outputs` (1e-10): ignored; profiles diverge ~1e-3.

## Related

- [Full forward-backward reconstruction proposal](../proposals/mugration-full-reconstruction-per-iteration.md)

## Related tickets

- [kb/tickets/mugration-iterative-gtr-golden-master-divergence.md](../tickets/mugration-iterative-gtr-golden-master-divergence.md)
