# Localize the residual ~1e-3 mugration marginal-confidence divergence from v0

With the GTR rate-optimizer fixed (`BrentBracketed`) and D1/D2 made opt-in (default v0), v1 reproduces v0 trait assignments for zika, zika-with-weights, and lassa. The remaining work is a single open problem: v1's marginal confidence profiles still differ from v0 by ~1e-3, even for unweighted, informative-root datasets where D1 and D2 are no-ops.

## Task

- Localize the ~1e-3 marginal-confidence divergence (example: zika_20_country root `NODE_0000000`, v0 `0.4812` vs v1 `0.4800`). Candidates: GTR `W`/eigendecomposition, branch-length-to-GTR mapping, normalization order in the marginal passes.
- Drive the divergence below the 1e-6 oracle tolerance, then un-ignore the affected cases:
  - `test_gm_mugration_outputs_v1_divergence` (dengue, tb, rsv, mpox) -> fold into `test_gm_mugration_outputs`.
  - `test_gm_mugration_confidence_zika` (1e-6) and `test_gm_mugration_confidence_outputs` (1e-10).
- NEVER widen these tolerances to pass; the gap is a real numerical defect.

## Related issues

- Source: [kb/issues/M-mugration-iterative-gtr.md](../issues/M-mugration-iterative-gtr.md)
