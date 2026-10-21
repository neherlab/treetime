# Mugration golden master parity with v0

v1 refines the mugration GTR model in the order of v0's `reconstruct_discrete_traits()` ([packages/legacy/treetime/treetime/wrappers.py#L785-L811](../../packages/legacy/treetime/treetime/wrappers.py#L785-L811)). With this order, v1 reproduces the v0 trait assignments of all seven golden-master datasets. The confidence profiles agree with v0 to at most 2.5e-8, but four datasets still exceed the 1e-10 tolerance of the confidence test.

## Refinement order

`fn refine_gtr()` [packages/treetime/src/gtr/refinement.rs](../../packages/treetime/src/gtr/refinement.rs) with the rate-optimizing policy runs the same steps as v0:

1. A marginal update under the start model, and the first fit from it
2. A full marginal update under the first fitted model (v0 `infer_ancestral_sequences(infer_gtr=True)`, [packages/legacy/treetime/treetime/treeanc.py#L564-L566](../../packages/legacy/treetime/treetime/treeanc.py#L564-L566)), then the first rate step
3. For each of the `iterations`: a fit, then a rate step. Each fit counts with the forward messages and the root state of the update of step 2, and with the backward messages of the last rate step. In v0, `infer_gtr()` reuses the existing marginal reconstruction ([packages/legacy/treetime/treetime/treeanc.py#L1544](../../packages/legacy/treetime/treetime/treeanc.py#L1544)), and `optimize_gtr_rate()` runs only backward passes and does not assign `root.cseq` ([packages/legacy/treetime/treetime/treeanc.py#L1685-L1689](../../packages/legacy/treetime/treetime/treeanc.py#L1685-L1689)), so the root state count ([packages/legacy/treetime/treetime/treeanc.py#L1605-L1613](../../packages/legacy/treetime/treetime/treeanc.py#L1605-L1613)) reads the root state of step 2
4. The sampling-bias scale and a final marginal update

## Resolved contributors

- **Forward messages and root state of the later fits**: v1 counted every fit with the forward messages of the start-model update, and the root state with the posterior of the last rate step. Both now come from the full update under the first fitted model, as in v0. This change moved every confidence profile toward v0 by 4 to 12 orders of magnitude (see Measurements)
- **Node names in the golden-master test**: v0 ladderizes the tree before it names internal nodes ([packages/legacy/treetime/treetime/treeanc.py#L446-L480](../../packages/legacy/treetime/treetime/treeanc.py#L446-L480)), and v1 names them in file order. On the trees of dengue, tb, rsv, and mpox, the same `NODE_*` name pointed to different nodes, so the earlier failures of these cases (trait mismatches and confidence differences near 1) compared different nodes. The test now gives the v1 nodes the names that v0 assigns
- **GTR rate optimizer**: the rate optimizer briefly used argmin `BrentOpt` (golden-section seeded), which converged to a different rate than v0's interior-seeded scipy `brent` and flipped assignments at ambiguous nodes. Fixed by `BrentBracketed` at [packages/treetime/src/gtr/brent_bracketed.rs](../../packages/treetime/src/gtr/brent_bracketed.rs)
- **Initial-pi pseudo-count and uninformative-root filtering**: these are opt-in flags defaulting to v0 (`--smooth-initial-pi`, `--filter-uninformative-root`). They do not change the default reconstruction

## Measurements

Command: `./dev/docker/run just test-rs -E 'test(/test_gm_mugration/)' --run-ignored all`. The columns give the number of nodes whose trait differs from v0 and the largest absolute difference of a confidence value over all nodes, with nodes matched by the names that v0 assigns. "Before" is the start-model order, "after" is the v0 order.

| Case                       | Before: traits | Before: confidence | After: traits | After: confidence | Confidence test (1e-10) |
| -------------------------- | -------------- | ------------------ | ------------- | ----------------- | ----------------------- |
| `zika_20_country`          | 0              | 1.3e-2             | 0             | 3.9e-14           | pass (also at 1e-6)     |
| `zika_20_country_weights`  | 0              | 3.1e-3             | 0             | 2.5e-8            | fail                    |
| `dengue_20_country`        | 1              | 6.5e-2             | 0             | 1.1e-13           | pass                    |
| `tb_20_cluster`            | 0              | 2.5e-3             | 0             | 7.4e-11           | pass                    |
| `lassa_L_20_country`       | 0              | 2.9e-3             | 0             | 5.8e-10           | fail                    |
| `rsv_a_20_country`         | 0              | 3.2e-2             | 0             | 3.0e-10           | fail                    |
| `mpox_clade_ii_20_country` | 0              | 1.6e-2             | 0             | 5.8e-9            | fail                    |

Before both changes, 3 of 14 test cases passed: the trait assignments of zika, zika with weights, and lassa. After both changes, 10 of 14 pass.

## Remaining divergence

The confidence profiles of `zika_20_country_weights`, `lassa_L_20_country`, `rsv_a_20_country`, and `mpox_clade_ii_20_country` differ from v0 by 3.0e-10 to 2.5e-8. Per project rules, a numerical error above 1e-6 against the v0 oracle is a defect; these differences are below that, but above the 1e-10 tolerance of `test_gm_mugration_confidence_outputs_v1_divergence`. The origin is not localized.

Candidate causes, not yet tested against the golden masters:

- **Evaluation point of the rate step**: v0's `optimize_gtr_rate()` sets the rate to the optimum without a new backward pass, so the backward messages that its next fit counts with come from the last Brent evaluation; v1 recomputes them at the optimum
- **Rate precision of the Brent search**: the stop rule and tolerance of `BrentBracketed` against scipy `brent`
- **Rate restore without a bracket**: `RateRangeEdge::KeepStart` ([M-mugration-gtr-rate-restore-inconsistency.md](M-mugration-gtr-rate-restore-inconsistency.md))
- **Argmax ties**: a tied root state follows rounding ([N-mugration-tied-state-follows-rounding.md](N-mugration-tied-state-follows-rounding.md)); this changes the root state count only at a tie
- **GTR numerics**: the exchangeability matrix `W` and its eigendecomposition; the branch length passed to `expQt`; the normalization order in the marginal passes

## Leaf evidence

Discrete leaves keep their observed trait as evidence in every pass, as v0 does, so a leaf posterior never feeds back as evidence. The golden-master datasets have no leaf with a missing trait, so their outputs do not depend on this rule.

> [!WARNING]
> **Needs review.** Measured before the v0 refinement order: on `data/zika/20` with the country of three leaves set to `?`, v1 fitted the rate 7.7146 and v0 fitted 6.4068. Measure again on the current code.

## Affected golden master tests

[packages/treetime/src/mugration/__tests__/test_gm_mugration.rs](../../packages/treetime/src/mugration/__tests__/test_gm_mugration.rs):

- `test_gm_mugration_outputs`: all seven datasets pass
- `test_gm_mugration_confidence_zika` (1e-6): passes
- `test_gm_mugration_confidence_outputs` (1e-10): dengue and tb pass
- `test_gm_mugration_confidence_outputs_v1_divergence` (1e-10): zika with weights, lassa, rsv, and mpox are ignored

## Fix

- Localize the remaining difference, starting from the largest one, `zika_20_country_weights` (2.5e-8), and test the candidate causes above
- Drive the difference below 1e-10, then fold the cases of `test_gm_mugration_confidence_outputs_v1_divergence` into `test_gm_mugration_confidence_outputs`
- Keep the test tolerances unchanged: the gap is a numerical defect, so a wider tolerance would hide it

## Related

- [kb/proposals/mugration-full-reconstruction-per-iteration.md](../proposals/mugration-full-reconstruction-per-iteration.md)
- [M-mugration-gtr-rate-restore-inconsistency.md](M-mugration-gtr-rate-restore-inconsistency.md)
- [N-mugration-tied-state-follows-rounding.md](N-mugration-tied-state-follows-rounding.md)
