# Timetree validation data

Simulated trees and sequences with known node dates, used to validate TreeTime timetree inference.

`./generate_all.sh [n_reps=10] [n_jobs=available CPUs]` writes `data/`. One replicate: `python simulate.py --coal exp --rep 0`
(with `PYTHONPATH=../packages/legacy/treetime`).

## Trees (`coalescent.py`)

Kingman coalescent, simulated backwards in time. `tau` is years before the most recent sample; a pair of lineages
merges at rate `1/Tc(tau)`.

| schedule    | Tc(tau) in years                              |
| ----------- | --------------------------------------------- |
| `const`     | 10                                            |
| `exp`       | 10 exp(-tau/5)                                |
| `piecewise` | 10 for tau < 5, 1 for 5 <= tau < 6, 3 beyond  |

Sampling: 10 tips at each of tau = 0, 0.1, ..., 3.9 (400 tips). The most recent sample is dated 2025.0.

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

## Validation runs (`run.py`, `evaluate.py`)

```
python run.py --suite ancestral --set core      # also: --set smoke (rep00, L=1000)
python run.py --suite timetree --set core
python run.py --suite timetree --set extra      # rep00-02, L=1000, mu in {0.003, 0.001}
python evaluate.py --suite timetree             # writes results/timetree.tsv, prints a summary
```

Jobs run in parallel (`-j`, default: available CPUs), skip finished runs, and record command, return code, run time
and peak memory in `run.json`. `--binary` selects the treetime executable (default `../target/release/treetime`),
`--filter` restricts jobs by a regex on the output path.

Inputs (`inputs.py`, written to `data/.../inputs/`) are divergence trees on the true topology: `exact` (mu t),
`realized` (true parent-child differences / L) and `rate` (realized times lognormal(0, 0.25) per branch).

- `ancestral`: `exact` tree; marginal, marginal with `--gtr-iterations 3`, parsimony. Metrics: per-site error of
  internal sequences (by node age), error of pi, W and mu of the inferred GTR.
- `timetree`: `realized` tree, alignment, `--clock-filter 0`; core variants `<root>-<prior>-<clock>` with root
  `keep`/`reroot`, coalescent prior `none`/`opt`/`skyline`, clock rate `est`/`fixed` (true mu). Metrics: clock rate,
  root date, internal-node date bias and RMSE (by node age), root recovery, Tc error and coverage, 90% CI coverage
  (`extra` set), and the number of nodes dated before their parent.

Nodes are matched to the truth by their set of descendant tips.
