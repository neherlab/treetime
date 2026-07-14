# Mugration golden master parity with v0

v1 implements iterative GTR inference for mugration, matching v0's `reconstruct_discrete_traits()` ([packages/legacy/treetime/treetime/wrappers.py#L653-L811](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811)) algorithm structure. With the default (v0) inference policy, v1 reproduces v0 trait assignments for zika, zika-with-weights, and lassa. The remaining datasets (dengue, tb, rsv, mpox) still diverge at a few ambiguous internal nodes, and all confidence profiles diverge from v0 by ~1e-3.

## Resolved contributors

- GTR rate-optimizer regression: the rate optimizer briefly used argmin `BrentOpt` (golden-section seeded), which converged to a different rate than v0's interior-seeded scipy `brent` and flipped assignments at ambiguous nodes. Fixed by `BrentBracketed` at [packages/treetime/src/gtr/brent_bracketed.rs](../../packages/treetime/src/gtr/brent_bracketed.rs).
- D1 (initial-pi pseudo-count) and D2 (uninformative-root filtering) are now opt-in flags defaulting to v0 (`--smooth-initial-pi`, `--filter-uninformative-root`). They no longer perturb the default reconstruction.

## Remaining divergence: residual ~1e-3 marginal-confidence difference

The unresolved cause is a ~1e-3 difference in v1's marginal confidence profiles relative to v0, present even for unweighted, informative-root datasets where D1 and D2 are no-ops. Example (zika_20_country root `NODE_0000000`): v0 `0.4812`, v1 `0.4800`. This small difference is below the old ~1.34e-2 figure but still tips the argmax at near-tied nodes, which is why dengue/tb/rsv/mpox assignments differ. The origin (marginal/GTR numerics) is not yet localized.

Per project rules, numerical error > 1e-6 against the v0 oracle is a defect. This remains open.

## Affected golden master tests

- `test_gm_mugration_outputs`: zika, zika_weights, lassa pass; dengue, tb, rsv, mpox ignored (`test_gm_mugration_outputs_v1_divergence`).
- `test_gm_mugration_confidence_zika` (1e-6) and `test_gm_mugration_confidence_outputs` (1e-10): ignored; profiles diverge ~1e-3.

## Related

- [Full forward-backward reconstruction proposal](../proposals/mugration-full-reconstruction-per-iteration.md)

## Related tickets

- [kb/tickets/mugration-iterative-gtr-golden-master-divergence.md](../tickets/mugration-iterative-gtr-golden-master-divergence.md)
