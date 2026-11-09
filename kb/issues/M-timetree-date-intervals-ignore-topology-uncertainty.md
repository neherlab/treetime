# The date intervals ignore the uncertainty of the topology

> [!NOTE]
> **Deferred.** No direct fix exists without sampling alternative trees or resolutions of the polytomies. Left for later. A cheap heuristic to test first is listed under Fix direction.

## Problem

`--confidence` writes 90% date intervals that condition on the input topology ([packages/treetime/src/timetree/confidence.rs](../../packages/treetime/src/timetree/confidence.rs)). On a reconstructed topology at low temporal signal, the error of the topology is larger than the error of the dates, and the intervals do not cover the true date of a clade at the node that joins its tips.

Validation, FastTree input (`run.py --set recon`), skyline prior, estimated clock rate, the 90% interval of the inferred common ancestor of the tips of each true internal node ([kb/reports/2026-11-08-timetree-short-branch-collapse-on-true-and-reconstructed-topologies.md](../reports/2026-11-08-timetree-short-branch-collapse-on-true-and-reconstructed-topologies.md#date-intervals)), coverage clean / messy:

| $L\mu$ | 400 tips | 1200 tips | 4000 tips |
| --- | --- | --- | --- |
| 0.3 | 0.08 / 0.09 | 0.04 / 0.04 | 0.03 / 0.02 |
| 1 | 0.30 / 0.35 | 0.15 / 0.18 | 0.09 / 0.10 |
| 3 | 0.56 / 0.57 | 0.46 / 0.45 | 0.33 / 0.32 |
| 10 | 0.76 / 0.75 | 0.69 / 0.69 | 0.59 / 0.58 |
| 30 | 0.86 / 0.85 | 0.82 / 0.81 | 0.76 / 0.76 |
| 100 | | 0.85 / 0.85 | 0.85 / 0.85 |

The values are for the default run with the short-branch collapse; without it they differ by at most 0.12, lower at $L\mu \le 3$ on messy data and slightly higher at $L\mu \ge 10$. The coverage falls with the size of the tree, because a larger tree has more wrong or unresolved splits at a given $L\mu$. On the true topology the same intervals cover 0.90 without the collapse.

## Fix direction

- Proper: a sensitivity band over a set of input trees (bootstrap or near-optimal trees) or over random resolutions of the polytomies, as sketched in [kb/reports/stochastic.simulation.md](../reports/stochastic.simulation.md#66-topology-uncertainty)
- Heuristic, untested: widen the interval of a node by a factor $(n_\text{children}+1)/3$, which is 1 for a binary node and grows with the size of a polytomy. Test against `mrca_ci90_coverage` on the `recon` set and against `ci90_coverage` on the true topology, where it must not over-cover
- At least: state in the documentation of `--confidence` that the intervals condition on the input topology

## Validation

`mrca_ci90_coverage` of the `ft-skyline-ci-*` variants of the `recon` set, per $L\mu$ and tree size.
