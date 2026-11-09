# Timetree short-branch collapse on true and reconstructed topologies

Date: 2026-11-08

The pre-loop short-branch collapse ([kb/features/timetree.md](../features/timetree.md), `fn collapse_short_internal_edges()` in [packages/treetime/src/optimize/topology/collapse.rs](../../packages/treetime/src/optimize/topology/collapse.rs)) helps at low and medium temporal signal when the input topology is reconstructed, and costs accuracy and interval calibration when it is the true one. The validation suite dates only the true topology, so it measures the cost and not the benefit. Keep the collapse on by default; `--no-collapse-short-branches` turns it off.

## Setup

- Code: `claude-work` at `2a46c49a1`, which adds `--no-collapse-short-branches` to `timetree` (the setting `TimetreeParams::collapse_short_branches` existed; the CLI set it to `true`)
- Baseline: `f7565a6df`, before the pre-loop collapse
- Data: the validation sets, 400 tips (`data`, `data_messy`), 1200 tips (`data_n30`, `data_n30_messy`) and 4000 tips (`data_n100`, `data_n100_messy`); `const`, `exp` and `piecewise` coalescents, $L \in \{10^3, 10^4\}$, $\mu$ from $3\cdot10^{-4}$ to $10^{-2}$. Signal is given as $L\mu$, the expected substitutions per sequence per year
- True topology: the `realized` input (true substitutions per edge over $L$), which every timetree run of `run.py --set core|extra` uses
- Reconstructed topology: `run.py --set recon`, the `fasttree` input of [validation/inputs.py](../../validation/inputs.py) (`FastTree -nt -gtr` on the tip alignment that TreeTime reads), rerooted by least squares, with an estimated clock rate, each prior run with and without the collapse
- Scores, from [validation/evaluate.py](../../validation/evaluate.py):
  - matched: inferred internal nodes whose tip set is a true clade (`date_rmse`, `ci90_coverage`, `frac_nodes_matched`). It scores a different node set in each arm, and on a collapsed tree only the clades that survive
  - common ancestor: every true internal node scored at the inferred common ancestor of its tips (`mrca_rmse`, `mrca_ci90_coverage`). Both arms are scored on the same nodes; a collapsed clade is scored at its polytomy, and a clade broken by wrong splits at the deeper node that joins its tips
- Values are medians over sets; "paired" is the fraction of sets where the collapse gives the lower `mrca_rmse`
- Only runs with a tree prior (`--coalescent opt` or `skyline`) count as evidence. Without a prior the dates are not well determined, so a smaller error there means little
- $T_c$ error is `tc_log_mae`, the mean absolute log ratio to the true $T_c$ over the coalescent segments. For `opt`, which fits one constant $T_c$, it is reported on the `const` sets only: a constant has no true value when the true $T_c$ changes in time (`exp`, `piecewise`)

## True topology

The collapse removes true edges without substitutions, the edges of $0/L$ length. With it, 13% of the inferred internal nodes are true clades at $L\mu=0.3$, 37% at 1, 69% at 3, 89% at 10 (400 tips). Without it the dates match the baseline to the third digit at every prior and $L\mu$, so the collapse explains the whole timetree difference to `f7565a6df`.

400 tips, clean data, collapse / no collapse:

| prior | $L\mu$ | matched rmse (y) | $T_c$ log MAE | ci90 coverage |
| --- | --- | --- | --- | --- |
| opt | 0.3 | 0.92 / 0.34 | 1.92 / 0.14 (`const`) | |
| opt | 1 | 0.36 / 0.19 | 0.84 / 0.06 (`const`) | |
| opt | 3 | 0.19 / 0.14 | 0.16 / 0.06 (`const`) | |
| skyline | 0.3 | 1.04 / 0.24 | 2.26 / 0.38 | |
| skyline | 1 | 0.39 / 0.18 | 0.98 / 0.32 | |
| skyline, `--confidence`, $\mu=10^{-3}$ | 1 | | | 0.745 / 0.900 |
| skyline, `--confidence`, $\mu=3\cdot10^{-3}$ | 3 | | | 0.793 / 0.903 |

- At 400 tips the arms agree from $L\mu=10$ up
- Without a prior the collapse gives a smaller root error (6.3 against 15.9 y at $L\mu=0.3$), which does not count as evidence (see Setup)
- With the collapse the intervals are wider on average (0.50 against 0.35 y at $\mu=10^{-3}$) and still cover less: the dates of the polytomies are biased, not only too narrow
- Messy data shows the same pattern. Larger trees show a larger cost that reaches higher signal, because they have more true edges without a substitution: at 4000 tips 4% of the inferred internal nodes are true clades at $L\mu=0.3$, 34% at 3, 72% at 10

| tips | ci90 coverage, $\mu=10^{-3}$ | opt $T_c$ log MAE (`const`), $L\mu$ = 1 / 3 / 10 |
| --- | --- | --- |
| 400 | 0.745 / 0.900 | 0.84 / 0.06, 0.16 / 0.06, 0.05 / 0.05 |
| 1200 | 0.69 / 0.90 | 1.48 / 0.09, 0.49 / 0.05, 0.04 / 0.04 |
| 4000 | 0.62 / 0.88 | 2.41 / 0.10, 1.07 / 0.03, 0.14 / 0.02 |

- Ancestral and optimize results are identical between the code and the baseline at every tree size

## Reconstructed topology

FastTree resolves every polytomy, mostly by edges of near-zero length. On `const/rep00`, $L=10^3$, $\mu=10^{-3}$, it returns 245 splits, 108 of them true, and 142 of length below $10^{-5}$. The collapse removes most of the wrong splits: at $L\mu=3$, 68% of the surviving internal nodes are true clades against 52% without it, and 89% against 76% at $L\mu=10$.

400 tips, clean data, collapse / no collapse:

| prior | $L\mu$ | common-ancestor rmse (y) | paired | abs. root error (y) | $T_c$ log MAE |
| --- | --- | --- | --- | --- | --- |
| opt | 0.3 | 1.78 / 3.66 | 0.90 | 6.5 / 13.6 | 1.82 / 1.77 (`const`) |
| opt | 1 | 0.76 / 1.04 | 1.00 | 3.8 / 5.9 | 0.83 / 0.99 (`const`) |
| opt | 3 | 0.25 / 0.36 | 0.98 | 1.0 / 2.1 | 0.14 / 0.34 (`const`) |
| skyline | 0.3 | 1.87 / 2.01 | 0.73 | 3.8 / 4.1 | 2.28 / 2.00 |
| skyline | 1 | 0.78 / 0.88 | 0.93 | 1.4 / 1.8 | 0.97 / 0.81 |
| skyline | 3 | 0.26 / 0.32 | 0.97 | 1.1 / 1.1 | 0.42 / 0.45 |
| opt | 30 | 0.053 / 0.052 | 0.13 | 0.23 / 0.20 | 0.05 / 0.05 (`const`) |
| skyline | 30 | 0.054 / 0.053 | 0.13 | 0.28 / 0.20 | 0.32 / 0.32 |

- At $L\mu \le 3$ the collapse wins in 73–100% of the sets. Without it the dates of the clades are biased old (opt, $L\mu=1$: bias −0.65 against −0.40 y): the wrong splits get dates, and the common ancestor of a broken clade lies deeper
- At $L\mu=10$ the arms are about equal
- At $L\mu \ge 30$ the collapse loses most paired comparisons, by about 0.001 y in rmse, and the skyline root error grows from 0.20 to 0.28 y. Here most short FastTree edges are true, and $0.1/L$ removes some of them
- The skyline $T_c$ is slightly worse with the collapse at $L\mu \le 1$, and at $L\mu=0.3$ both arms miss $T_c$ by a factor of about 6 to 10
- Messy data shows larger gains: with the opt prior the collapse wins every set at $L\mu$ from 0.3 to 3 (rmse 1.94 against 3.27 y at $L\mu=0.3$, 0.70 against 1.12 at 1)
- Without a prior the collapse also gives smaller errors, which do not count as evidence (see Setup)
- 2 of 3840 runs fail, the same set (`piecewise/rep09`, $L=10^3$, $\mu=3\cdot10^{-4}$, messy, no prior) in both arms, at the grid point limit

### Larger trees

The same runs on 1200 and 4000 tips, prior variants only (`run.py --set recon --filter '^(?!.*/ft-none-)'`). Opt prior, collapse / no collapse:

| tips | data | $L\mu$ | common-ancestor rmse (y) | paired | abs. root error (y) | $T_c$ log MAE (`const`) |
| --- | --- | --- | --- | --- | --- | --- |
| 1200 | clean | 1 | 0.98 / 1.47 | 1.00 | 4.5 / 8.3 | 1.35 / 1.51 |
| 1200 | clean | 3 | 0.27 / 0.38 | 1.00 | 2.2 / 3.5 | 0.49 / 0.72 |
| 1200 | clean | 10 | 0.078 / 0.101 | 0.98 | 0.40 / 0.86 | 0.09 / 0.16 |
| 1200 | clean | 30 | 0.041 / 0.040 | 0.33 | 0.42 / 0.26 | 0.05 / 0.03 |
| 4000 | clean | 1 | 0.91 / 1.39 | 1.00 | 2.6 / 6.3 | 2.45 / 2.53 |
| 4000 | clean | 3 | 0.32 / 0.40 | 1.00 | 2.9 / 3.9 | 1.04 / 1.27 |
| 4000 | clean | 10 | 0.070 / 0.107 | 1.00 | 0.53 / 1.6 | 0.08 / 0.28 |
| 4000 | clean | 30 | 0.030 / 0.031 | 0.80 | 0.15 / 0.31 | 0.06 / 0.04 |
| 4000 | clean | 100 | 0.0156 / 0.0153 | 0.23 | 0.21 / 0.17 | 0.07 / 0.07 |
| 4000 | messy | 10 | 0.073 / 0.112 | 1.00 | 0.47 / 1.4 | 0.06 / 0.30 |

Skyline prior, collapse / no collapse, paired in parentheses:

| tips | data | $L\mu=0.3$ | $L\mu=1$ | $L\mu=3$ | $L\mu=10$ | $L\mu=30$ |
| --- | --- | --- | --- | --- | --- | --- |
| 1200 | clean | 2.20 / 2.39 (0.63) | 1.09 / 1.10 (0.47) | 0.28 / 0.35 (1.00) | 0.080 / 0.097 (0.97) | 0.042 / 0.040 (0.20) |
| 1200 | messy | 2.64 / 3.91 (0.80) | 1.23 / 1.29 (0.77) | 0.28 / 0.40 (1.00) | 0.081 / 0.098 (0.95) | 0.043 / 0.041 (0.23) |
| 4000 | clean | 2.15 / 2.74 (0.83) | 1.19 / 1.20 (0.43) | 0.34 / 0.39 (0.97) | 0.072 / 0.101 (1.00) | 0.030 / 0.032 (0.70) |
| 4000 | messy | 2.62 / 6.21 (1.00) | 1.28 / 1.34 (0.63) | 0.35 / 0.45 (0.98) | 0.074 / 0.108 (1.00) | 0.031 / 0.032 (0.57) |

- The range of signal where the collapse helps grows with the tree: up to $L\mu=3$ at 400 tips, up to 10 at 1200, and at 4000 tips it still wins most sets at 30. A larger tree has more short FastTree edges at a given $L\mu$, and more of them are wrong
- Where the collapse loses ($L\mu \ge 30$ at 1200 tips, 100 at 4000) the rmse differs by 0.002 y or less
- The collapse keeps more true clades among the inferred nodes: 0.71 against 0.54 at 4000 tips, $L\mu=10$, clean
- With the skyline prior at $L\mu=1$ the arms are about equal (43–77% paired), and the skyline $T_c$ is slightly worse with the collapse at $L\mu \le 1$, as at 400 tips
- 1 of 5760 runs fails: `const/rep02`, $L=10^3$, $\mu=10^{-3}$, messy, 1200 tips, skyline with `--confidence` and the collapse

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
- The coverage at the common ancestors falls with the size of the tree, in both arms: at $L\mu=10$, 0.76 at 400 tips, 0.69 at 1200, 0.59 at 4000 (collapse, clean). On the matched nodes at $L\mu=10$ the collapse gives 0.78 / 0.88 at 1200 tips and 0.75 / 0.87 at 4000

## Conclusions

- The collapse removes the arbitrary resolutions of a tree builder, and its benefit on reconstructed trees outweighs its cost on true ones. The benefit grows with the size of the tree and reaches $L\mu=30$ at 4000 tips. The default stays on
- The validation core and extra sets date the true topology, so they show only the cost. A result on them that depends on the collapse needs the `recon` set as well
- Open: the calibration of the dates of collapsed polytomies, and a threshold that keeps more true short edges at high signal

## Reproduction

From `validation/`, with the FastTree binary at `../.out/FastTree` (or named by `FASTTREE`):

```
python run.py --suite timetree --set recon --data data --results results_recon
python evaluate.py --suite timetree --data data --results results_recon
```

The larger trees used `--data data_n30|data_n100` (and `_messy`) with `--filter '^(?!.*/ft-none-)'`, which skips the runs without a prior.

The true-topology runs without the collapse used a wrapper that adds `--no-collapse-short-branches` to every `timetree` call, with `run.py --suite timetree --set core|extra --binary <wrapper>`.
