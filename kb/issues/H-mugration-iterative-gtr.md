# Mugration golden master parity with v0

v1 implements iterative GTR inference for mugration, matching v0's `reconstruct_discrete_traits()` ([packages/legacy/treetime/treetime/wrappers.py#L653-L811](../../packages/legacy/treetime/treetime/wrappers.py#L653-L811)) algorithm structure. With the default (v0) inference policy, v1 reproduces v0 trait assignments for zika, zika-with-weights, and lassa. The remaining datasets (dengue, tb, rsv, mpox) still diverge at a few ambiguous internal nodes, and the confidence profiles of every dataset diverge from v0, by up to 3.5e-2.

## Resolved contributors

- GTR rate-optimizer regression: the rate optimizer briefly used argmin `BrentOpt` (golden-section seeded), which converged to a different rate than v0's interior-seeded scipy `brent` and flipped assignments at ambiguous nodes. Fixed by `BrentBracketed` at [packages/treetime/src/gtr/brent_bracketed.rs](../../packages/treetime/src/gtr/brent_bracketed.rs).
- D1 (initial-pi pseudo-count) and D2 (uninformative-root filtering) are now opt-in flags defaulting to v0 (`--smooth-initial-pi`, `--filter-uninformative-root`). They no longer perturb the default reconstruction.

## Remaining divergence: marginal-confidence difference

v1's marginal confidence profiles differ from v0, also for unweighted, informative-root datasets where D1 and D2 are no-ops. A run of the ignored golden-master tests on `rust` at `3825359a2` (`--run-ignored only`) fails in all 11 cases:

- zika_20_country root `NODE_0000000`: v0 `0.4812`, v1 `0.4800`
- dengue_20_country: up to 3.5e-2 at one node (v0 `0.938`, v1 `0.903`)
- mpox_clade_ii_20_country: v0 `0.068`, v1 `0.073` in one profile
- lassa_L_20_country: differences of about 1e-8 in a profile near 1

The difference tips the argmax at near-tied nodes, so the trait assignments of dengue, tb, rsv, and mpox differ from v0. The origin is not yet localized.

Per project rules, numerical error > 1e-6 against the v0 oracle is a defect. This remains open.

## Candidate causes (not yet tested against the golden masters)

- **Stale forward messages**: v1 counts transitions in every iteration with the forward messages of the initial uniform model, while v0 counts with forward messages of the first fitted GTR. This is a verified divergence with its own issue: [M-mugration-gtr-refinement-counts-with-forward-messages-of-the-initial-model.md](M-mugration-gtr-refinement-counts-with-forward-messages-of-the-initial-model.md)
- **Rate-optimizer evaluation point**: v0's `optimize_gtr_rate` sets the rate to the optimum without a new backward pass, so its subtree messages come from the last Brent evaluation; v1 recomputes them at the optimum
- **GTR numerics**: the exchangeability matrix `W` and its eigendecomposition
- **Branch-length-to-GTR mapping**: the branch length the marginal passes and the transition counts pass to `expQt`
- **Normalization order** in the marginal passes

## Leaf evidence

Discrete leaves keep their observed trait as evidence in every pass ([packages/treetime/src/partition/marginal/discrete/partition.rs](../../packages/treetime/src/partition/marginal/discrete/partition.rs)), as v0 does, so a leaf posterior never feeds back as evidence. This is not the cause of the divergence above:

- The golden master datasets have no leaf with a missing trait. Their outputs do not depend on this rule, and zika_20_country root `NODE_0000000` is `0.4800` (v0 `0.4812`)
- Inputs with missing traits diverge more. On `data/zika/20` with the country of three leaves set to `?`, v1 fits the rate 7.7146 and v0 fits 6.4068. The cause is not localized

## Affected golden master tests

- `test_gm_mugration_outputs`: zika, zika_weights, lassa pass; dengue, tb, rsv, mpox ignored (`test_gm_mugration_outputs_v1_divergence`).
- `test_gm_mugration_confidence_zika` (1e-6) and `test_gm_mugration_confidence_outputs` (1e-10): ignored; profiles diverge ~1e-3.

## Fix

- Localize the ~1e-3 marginal-confidence divergence, starting from zika_20_country root `NODE_0000000` (v0 `0.4812`, v1 `0.4800`)
- Drive the divergence below the 1e-6 oracle tolerance, then enable the ignored cases: fold `test_gm_mugration_outputs_v1_divergence` (dengue, tb, rsv, mpox) into `test_gm_mugration_outputs`, and remove the `#[ignore]` from `test_gm_mugration_confidence_zika` (1e-6) and `test_gm_mugration_confidence_outputs` (1e-10) [packages/treetime/src/mugration/__tests__/test_gm_mugration.rs](../../packages/treetime/src/mugration/__tests__/test_gm_mugration.rs)
- Keep the test tolerances unchanged: the gap is a numerical defect, so a wider tolerance would hide it

## Related

- [kb/proposals/mugration-full-reconstruction-per-iteration.md](../proposals/mugration-full-reconstruction-per-iteration.md)
- [M-mugration-gtr-rate-restore-inconsistency.md](M-mugration-gtr-rate-restore-inconsistency.md)
- [M-mugration-gtr-refinement-counts-with-forward-messages-of-the-initial-model.md](M-mugration-gtr-refinement-counts-with-forward-messages-of-the-initial-model.md)
