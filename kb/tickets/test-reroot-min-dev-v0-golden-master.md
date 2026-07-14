# Add a v0 golden master for minimum-deviation rerooting

Capture the v0 minimum-deviation root on a small committed dataset and compare the complete v1 reroot result against that independent oracle.

## Acceptance criteria

- The capture records the selected root edge, split fraction, and resulting incident branch lengths from v0.
- The v1 optimize reroot path selects the same edge as v0 and matches the split fraction and resulting incident branch lengths with absolute error at most $10^{-6}$.
- The clock reroot comparison remains a diagnostic expected to expose the separately tracked objective and split-optimizer divergences; this test ticket does not change those behaviors.
- Endpoint roots and ties have explicit deterministic expectations.
- The expected values are captured from v0 and are not computed by v1 helpers.
- Perturbing the expected edge or split fraction makes the test fail clearly.

## Related issues

- Source: [kb/issues/N-reroot-missing-v0-golden-master.md](../issues/N-reroot-missing-v0-golden-master.md)
- [kb/issues/M-clock-mindev-wrong-objective.md](../issues/M-clock-mindev-wrong-objective.md)
- [kb/issues/N-reroot-split-optimizer-default-diverges-from-v0.md](../issues/N-reroot-split-optimizer-default-diverges-from-v0.md)
