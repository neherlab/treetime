# Timetree validation data

Simulated trees and sequences with known node dates, used to validate TreeTime timetree inference.

`./generate_all.sh [n_reps=10] [n_jobs=8]` writes `data/`. One replicate: `python simulate.py --coal exp --rep 0`
(with `PYTHONPATH=../packages/legacy/treetime`).

## Trees (`coalescent.py`)

Kingman coalescent, simulated backwards in time. `tau` is years before the most recent sample; a pair of lineages
merges at rate `1/Tc(tau)`.

| schedule    | Tc(tau) in years                              |
| ----------- | --------------------------------------------- |
| `const`     | 10                                            |
| `exp`       | 10 exp(-tau/5)                                |
| `piecewise` | 10 for tau < 5, 1 for 5 <= tau < 6, 3 beyond  |

Sampling: 100 tips at each of tau = 0, 0.1, ..., 3.9 (4000 tips). The most recent sample is dated 2025.0.

## Sequences (`simulate.py`)

Per (L, mu): random GTR on ACGT (pi uniform-random and normalised, W symmetric uniform-random) scaled so that `mu`
is substitutions per site and year; root sequence drawn from pi; evolved along the tree with legacy
`treetime.seqgen.SeqGen`. Each (coal, rep, L, mu) has its own deterministic seed.

`GTR.random` is not used: it defaults to the gap-containing alphabet, does not keep W symmetric, and `assign_rates`
multiplies `mu` by the raw average rate.

## Layout

```
data/<coal>/rep<NN>/
  tree.nwk            true time tree, named nodes, branch lengths in years
  tip_dates.tsv       name, date, tau   (TreeTime --dates format)
  node_dates.tsv      name, date, tau   for internal nodes incl. node_root
  L<L>_mu<mu>/
    tips.fasta.xz
    internal.fasta.xz
    gtr.json          mu, pi, W (normalised), L
```
