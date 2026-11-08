# The date intervals ignore the uncertainty of the topology

> [!IMPORTANT]
> **Decision required.** Whether `timetree` should report date intervals that include topology uncertainty, and how. Options: a sensitivity band over a set of input trees (bootstrap or near-optimal trees), as sketched in [kb/reports/stochastic.simulation.md](../reports/stochastic.simulation.md#66-topology-uncertainty); a statement in the documentation of `--confidence` that the intervals condition on the input topology; or both. The evidence below shows the size of the gap; it does not choose an option.

## Problem

`--confidence` writes 90% date intervals that condition on the input topology ([packages/treetime/src/timetree/confidence.rs](../../packages/treetime/src/timetree/confidence.rs)). On a reconstructed topology at low temporal signal, the error of the topology is larger than the error of the dates, and the intervals do not cover the true date of a clade at the node that joins its tips.

Validation, 400 tips, FastTree input (`run.py --set recon`), skyline prior, estimated clock rate, the 90% interval of the inferred common ancestor of the tips of each true internal node ([kb/reports/2026-11-08-timetree-short-branch-collapse-on-true-and-reconstructed-topologies.md](../reports/2026-11-08-timetree-short-branch-collapse-on-true-and-reconstructed-topologies.md#date-intervals)):

| $L\mu$ | coverage, clean | coverage, messy |
| --- | --- | --- |
| 0.3 | 0.08 | 0.09 |
| 1 | 0.30 | 0.35 |
| 3 | 0.56 | 0.57 |
| 10 | 0.76 | 0.75 |
| 30 | 0.86 | 0.85 |

The values are for the default run with the short-branch collapse; without it they differ by at most 0.08. On the true topology the same intervals cover 0.90 without the collapse.

## Fix direction

Open, see the decision above.

## Validation

`mrca_ci90_coverage` of the `ft-skyline-ci-*` variants of the `recon` set, per $L\mu$.
