# Timetree validation data

Simulated trees and sequences with known node dates, used to validate TreeTime timetree inference.

All scripts run on the host in the pixi environment of the repository (`pixi.toml`), from this directory.
`pixi run ./generate_all.sh [n_reps=10] [n_jobs=available CPUs] [tips="10 30 100"]` writes `data/` (10 tips per
sample time), `data_n30/` and `data_n100/`. One replicate:
`pixi run python simulate.py --coal exp --rep 0 [--tips-per-sample 30 --out data_n30]`. A 100-tip replicate takes about
3 minutes on one core.

## Trees (`coalescent.py`)

Kingman coalescent, simulated backwards in time. `tau` is years before the most recent sample; a pair of lineages
merges at rate `1/Tc(tau)`.

| schedule    | Tc(tau) in years                              |
| ----------- | --------------------------------------------- |
| `const`     | 10                                            |
| `exp`       | 10 exp(-tau/5)                                |
| `piecewise` | 10 for tau < 5, 1 for 5 <= tau < 6, 3 beyond  |

Sampling: 10, 30 or 100 tips at each of tau = 0, 0.1, ..., 3.9 (400, 1200 or 4000 tips). The most recent sample is
dated 2025.0. The seeds of the larger sets include the tip count, so their trees are independent of the 10-tip trees.

## Sequences (`simulate.py`)

Per (L, mu): random GTR on ACGT (pi uniform-random and normalised, W symmetric uniform-random) scaled so that `mu`
is substitutions per site and year; root sequence drawn from pi; evolved along the tree with legacy
`treetime.seqgen.SeqGen`. Each (coal, rep, L, mu) has its own deterministic seed.

`GTR.random` is not used: it defaults to the gap-containing alphabet, does not keep W symmetric, and `assign_rates`
multiplies `mu` by the raw average rate.

## Messy data (`degrade.py`)

`generate_all.sh` and `slurm.sh` also write a messy copy of every replicate: `data_messy/`, `data_n30_messy/`,
`data_n100_messy/`. One replicate: `pixi run python degrade.py --coal exp --rep 0 [--data data_n30]`. Each sequence
set gets

- gaps: 3 per 1000 sites, 3 bases long, each on a branch drawn in proportion to its length in time and inherited by all
  nodes below it (in `tips.fasta.xz` and `internal.fasta.xz`; events in `gaps.tsv`). Legacy SeqGen has no indels.
- ambiguous bases: each tip site (not under a gap) with probability 0.002 becomes a two- or three-base IUPAC code
  that contains the true base.
- runs of N: 5% of tips get one run of 5 to 100 Ns.

and 3% of tips lose their date, 5% keep only the month (`2021-10-XX`). Truth files keep their names; TreeTime gets
`tips_observed.fasta.xz` and `tip_dates_observed.tsv` (`run.py` uses `*_observed` files when present). Run, evaluate
and analyze it like a tip density: `--data data_messy --results results_messy`. Extra metrics: `gap_recall` and
`false_gap_rate` (ancestral; `site_error` counts only non-gap sites), and the date RMSE of tips with a missing or
month-only date (`tip_date_rmse_missing`, `tip_date_rmse_month`). Realized branch lengths do not count gaps.

The indel model of `optimize` is tested by pairing `fixed-topo` with `fixed-topo-noindel` on the same set
(`branches.tsv` has the gap events per branch, `n_indels`). With ρ gap events per substitution, s substitutions and k
gap events on a branch should become about (s + k) / (1 + ρ): branches with a gap grow by about one substitution, the
others shrink by 1 + ρ (ρ is 0.04-0.06 at mu = 0.0003). `analyze.py` reports this in `summary.md` and
`optimize_indel.png`. On clean data the two variants must give identical trees.

## Layout

```
data/<coal>/rep<NN>/        10 tips per sample time; data_n30/ and data_n100/ have the same layout
  tree.nwk            true time tree, named nodes, branch lengths in years
  tip_dates.tsv       name, date, tau   (TreeTime --metadata format)
  node_dates.tsv      name, date, tau   for internal nodes incl. node_root
  L<L>_mu<mu>/
    tips.fasta.xz
    internal.fasta.xz
    gtr.json          mu, pi, W (normalised), L
data_messy/...        as data/, plus tip_dates_observed.tsv, L*/tips_observed.fasta.xz and L*/gaps.tsv
```

## Validation runs (`run.py`, `evaluate.py`)

Build the binary first with `./dev/docker/run just build-release treetime`, which writes `../.out/treetime`. The runs
call it from the host, because the Python image has an older C library than the binary needs.

```
pixi run python run.py --suite ancestral --set core      # also: --set smoke (rep00, L=1000)
pixi run python run.py --suite timetree --set core
pixi run python run.py --suite timetree --set extra      # rep00-02, L=1000, mu in {0.003, 0.001}
pixi run python run.py --suite optimize --set core       # branch length optimization
pixi run python evaluate.py --suite timetree             # writes results/timetree.tsv, prints a summary
pixi run python analyze.py                               # results/report/: summary.md and figures
pixi run python analyze.py --baseline ../old/results     # also compare.md: per-run changes against an earlier run
pixi run python run.py --suite optimize --set core --data data_n100 --results results_n100   # larger trees: own data and results
pixi run python evaluate.py --suite optimize --data data_n100 --results results_n100
pixi run python analyze.py --results results_n100
```

On a Slurm cluster, `./slurm.sh` submits the whole pipeline (generate, run, evaluate, analyze for every tip density) as
dependent jobs, with run tasks of one sequence set each (`run.py --shard`). Build `../target/release/treetime` and
submit from the head node: the jobs run the Python of the conda environment `CONDA_ENV` (default `nextstrain`),
resolved at submission, and need no network. Resources and the partition are set by environment variables (see the
script header), e.g. `SBATCH_PARTITION=cpu RUN_CPUS=8 ./slurm.sh`. Logs go to `slurm_logs/`.

At 4000 tips (`data_n100`, L=10000, mu=0.01), one run takes about 10 s for `optimize` and 1.5 min for `timetree`, with
1.7 GB peak memory; set `-j` so that jobs times 2 GB fits in memory. `fixed-topo-dense` needs more than 8 GB there;
skip it with `--filter '^(?!.*fixed-topo-dense)'`.

Failed runs carry their TreeTime error message in the `error` column of the TSV.

Jobs run in parallel (`-j`, default: available CPUs), skip finished runs, and record command, return code, run time
and peak memory in `run.json`. `--binary` selects the treetime executable (default `../.out/treetime`),
`--filter` restricts jobs by a regex on the output path.

Inputs (`inputs.py`, written to `data/.../inputs/`) are divergence trees on the true topology: `exact` (mu t),
`realized` (true parent-child differences / L), `rate` (realized times lognormal(0, 0.25) per branch) and `perturbed`
(exact times lognormal(0, 0.3) per branch).

- `ancestral`: `exact` tree; marginal, parsimony. Metrics: per-site error of internal sequences (by node age), error
  of pi, W and mu of the inferred GTR.
- `timetree`: `realized` tree, alignment, `--clock-filter 0`; core variants `<root>-<prior>-<clock>` with root
  `keep`/`reroot`, coalescent prior `none`/`opt`/`skyline`, clock rate `est`/`fixed` (true mu). Metrics: clock rate,
  root date, internal-node date bias and RMSE (by node age), root recovery, Tc error and coverage, 90% CI coverage
  (`extra` set), and the number of nodes dated before their parent.

- `optimize`: `perturbed` tree, alignment; variants `fixed-topo` (no collapsing, merging or flipping of branches),
  `fixed-topo-dense`, `fixed-topo-noindel` (`--no-indels`) and `default`. Metrics, in substitutions (branch length times L): per-branch RMSE against the
  realized branch lengths (also for the perturbed input, and by realized substitutions), total length against realized
  and expected, zero and non-zero branches kept as such, squared deviation from mu t relative to Poisson noise, and the
  sum of the two root branches, which is the only identifiable quantity at the root. Per-branch lengths go to
  `branches.tsv` in each run directory; `optimize_dispersion.png` bins branches by mu L t and plots the variance of
  the optimized counts against the mean (Poisson: equal). Dense takes minutes per run;
  `--filter '^(?!.*fixed-topo-dense/.*(rep0[3-9]|L10000))'` limits it to rep00-02 at L=1000.

Nodes and branches are matched to the truth by their set of descendant tips.
