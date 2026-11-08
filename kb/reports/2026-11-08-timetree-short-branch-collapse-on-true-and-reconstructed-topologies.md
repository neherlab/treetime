# Timetree short-branch collapse on true and reconstructed topologies

Date: 2026-11-08

The pre-loop short-branch collapse ([kb/features/timetree.md](../features/timetree.md), `fn collapse_short_internal_edges()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs)) helps at low and medium temporal signal when the input topology is reconstructed, and costs accuracy and interval calibration when it is the true one. The validation suite dates only the true topology, so it measures the cost and not the benefit. Keep the collapse on by default; `--no-collapse-short-branches` turns it off.

## Setup

- Code: `claude-work` at `2a46c49a1`, which adds `--no-collapse-short-branches` to `timetree` (the setting `TimetreeParams::collapse_short_branches` existed; the CLI set it to `true`)
- Baseline: `f7565a6df`, before the pre-loop collapse
- Data: the validation sets, 400 tips (`data`, `data_messy`) and 1200 tips (`data_n30`, `data_n30_messy`); `const`, `exp` and `piecewise` coalescents, $L \in \{10^3, 10^4\}$, $\mu$ from $3\cdot10^{-4}$ to $10^{-2}$. Signal is given as $L\mu$, the expected substitutions per sequence per year
- True topology: the `realized` input (true substitutions per edge over $L$), which every timetree run of `run.py --set core|extra` uses
- Reconstructed topology: `run.py --set recon`, the `fasttree` input of [validation/inputs.py](../../validation/inputs.py) (`FastTree -nt -gtr` on the tip alignment that TreeTime reads), rerooted by least squares, with an estimated clock rate, each prior run with and without the collapse
- Scores, from [validation/evaluate.py](../../validation/evaluate.py):
  - matched: inferred internal nodes whose tip set is a true clade (`date_rmse`, `ci90_coverage`, `frac_nodes_matched`). It scores a different node set in each arm, and on a collapsed tree only the clades that survive
  - common ancestor: every true internal node scored at the inferred common ancestor of its tips (`mrca_rmse`, `mrca_ci90_coverage`). Both arms are scored on the same nodes; a collapsed clade is scored at its polytomy, and a clade broken by wrong splits at the deeper node that joins its tips
- Values are medians over sets; "paired" is the fraction of sets where the collapse gives the lower `mrca_rmse`

## True topology

The collapse removes true edges without substitutions, the edges of $0/L$ length. With it, 13% of the inferred internal nodes are true clades at $L\mu=0.3$, 37% at 1, 69% at 3, 89% at 10 (400 tips). Without it the dates match the baseline to the third digit at every prior and $L\mu$, so the collapse explains the whole timetree difference to `f7565a6df`.

400 tips, clean data, collapse / no collapse:

| prior | $L\mu$ | matched rmse (y) | $T_c$ abs. error | ci90 coverage |
| --- | --- | --- | --- | --- |
| opt | 0.3 | 0.92 / 0.34 | 2.00 / 0.22 | |
| opt | 1 | 0.36 / 0.19 | 0.98 / 0.14 | |
| opt | 3 | 0.19 / 0.14 | 0.29 / 0.10 | |
| skyline | 0.3 | 1.04 / 0.24 | 2.26 / 0.38 | |
| skyline | 1 | 0.39 / 0.18 | 0.98 / 0.32 | |
| skyline, `--confidence`, $\mu=10^{-3}$ | 1 | | | 0.745 / 0.900 |
| skyline, `--confidence`, $\mu=3\cdot10^{-3}$ | 3 | | | 0.793 / 0.903 |

- From $L\mu=10$ up the arms agree
- Without a prior the collapse gives a smaller root error (6.3 against 15.9 y at $L\mu=0.3$) and a smaller matched rmse, but the matched rmse covers only the 13% of nodes that survive
- With the collapse the intervals are wider on average (0.50 against 0.35 y at $\mu=10^{-3}$) and still cover less: the dates of the polytomies are biased, not only too narrow
- Messy data and 1200 tips show the same pattern. At 1200 tips, messy, ci90 coverage is 0.75 against 0.885 at $\mu=10^{-3}$, and the opt $T_c$ error 1.47 against 0.08 at $L\mu=1$
- Ancestral and optimize results are identical between the code and the baseline at both tree sizes

## Reconstructed topology

FastTree resolves every polytomy, mostly by edges of near-zero length. On `const/rep00`, $L=10^3$, $\mu=10^{-3}$, it returns 245 splits, 108 of them true, and 142 of length below $10^{-5}$. The collapse removes most of the wrong splits: at $L\mu=3$, 68% of the surviving internal nodes are true clades against 52% without it, and 89% against 76% at $L\mu=10$.

400 tips, clean data, collapse / no collapse:

| prior | $L\mu$ | common-ancestor rmse (y) | paired | abs. root error (y) | $T_c$ log MAE |
| --- | --- | --- | --- | --- | --- |
| none | 0.3 | 5.9 / 6.9 | 0.77 | 6.2 / 14.7 | |
| none | 1 | 1.32 / 2.13 | 1.00 | 2.7 / 5.4 | |
| none | 3 | 0.41 / 0.64 | 0.98 | 1.5 / 2.0 | |
| opt | 0.3 | 1.78 / 3.66 | 0.90 | 6.5 / 13.6 | 1.88 / 2.77 |
| opt | 1 | 0.76 / 1.04 | 1.00 | 3.8 / 5.9 | 0.89 / 1.17 |
| opt | 3 | 0.25 / 0.36 | 0.98 | 1.0 / 2.1 | 0.22 / 0.38 |
| skyline | 1 | 0.78 / 0.88 | 0.93 | 1.4 / 1.8 | 0.97 / 0.81 |
| skyline | 3 | 0.26 / 0.32 | 0.97 | 1.1 / 1.1 | 0.42 / 0.45 |
| opt | 30 | 0.053 / 0.052 | 0.13 | 0.23 / 0.20 | 0.17 / 0.16 |
| skyline | 30 | 0.054 / 0.053 | 0.13 | 0.28 / 0.20 | 0.32 / 0.32 |

- At $L\mu \le 3$ the collapse wins in 77–100% of the sets. Without it the dates of the clades are biased old (opt, $L\mu=1$: bias −0.65 against −0.40 y): the wrong splits get dates, and the common ancestor of a broken clade lies deeper
- At $L\mu=10$ the arms are about equal
- At $L\mu \ge 30$ the collapse loses most paired comparisons, by about 0.001 y in rmse, and the skyline root error grows from 0.20 to 0.28 y. Here most short FastTree edges are true, and $0.1/L$ removes some of them
- Messy data shows larger gains: without a prior at $L\mu=0.3$, rmse 5.5 against 10.9 y and root error 6.9 against 21.9 y; with the opt prior the collapse wins every set at $L\mu$ from 0.3 to 3
- 2 of 3840 runs fail, the same set (`piecewise/rep09`, $L=10^3$, $\mu=3\cdot10^{-4}$, messy, no prior) in both arms, at the grid point limit

## Date intervals

Skyline prior with `--confidence`, FastTree input, clean data, collapse / no collapse:

| $L\mu$ | ci90 coverage, matched nodes | ci90 coverage, common ancestors |
| --- | --- | --- |
| 0.3 | 0.52 / 0.57 | 0.08 / 0.09 |
| 1 | 0.73 / 0.80 | 0.30 / 0.33 |
| 3 | 0.76 / 0.87 | 0.56 / 0.56 |
| 10 | 0.81 / 0.87 | 0.76 / 0.78 |
| 30 | 0.87 / 0.89 | 0.86 / 0.88 |

- On the matched nodes the collapse lowers the coverage, as on the true topology: [kb/issues/M-timetree-short-branch-collapse-lowers-date-interval-coverage.md](../issues/M-timetree-short-branch-collapse-lowers-date-interval-coverage.md)
- At the common ancestors the coverage is far below 0.9 in both arms at low signal. The intervals condition on one topology, and at $L\mu \le 3$ the topology error is larger than the dating error: [kb/issues/M-timetree-date-intervals-ignore-topology-uncertainty.md](../issues/M-timetree-date-intervals-ignore-topology-uncertainty.md)

## Conclusions

- The collapse removes the arbitrary resolutions of a tree builder, and its benefit on reconstructed trees outweighs its cost on true ones. The default stays on
- The validation core and extra sets date the true topology, so they show only the cost. A result on them that depends on the collapse needs the `recon` set as well
- Open: the calibration of the dates of collapsed polytomies, and a threshold that keeps more true short edges at high signal

## Reproduction

From `validation/`, with the FastTree binary at `../.out/FastTree` (or named by `FASTTREE`):

```
python run.py --suite timetree --set recon --data data --results results_recon
python evaluate.py --suite timetree --data data --results results_recon
```

The true-topology runs without the collapse used a wrapper that adds `--no-collapse-short-branches` to every `timetree` call, with `run.py --suite timetree --set core|extra --binary <wrapper>`.
