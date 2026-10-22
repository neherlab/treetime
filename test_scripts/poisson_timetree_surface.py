"""Clock rate and node times of a fixed tree under Poisson branch likelihoods and the exact Kingman coalescent density.

Purpose
  A reference for the clock rate that maximizes the marginal likelihood Z(mu) of an exact model, against which the
  rate step of TreeTime can be compared. The `mcmc` command samples the node times and estimates the slope
  d log Z / d ln mu; the maximum of Z lies where the slope crosses zero.

Model
  log p(k, t | mu) = sum_b [k_b ln(mu L tau_b) - mu L tau_b - ln k_b!] + log p_Kingman(t | Tc)

  k_b are the substitutions on branch b (branch length times L, rounded), tau_b = t_child - t_parent in years, and
  L is the alignment length. The Kingman density of heterochronous samples is
  log p_Kingman = -M ln Tc - int k(s) (k(s) - 1) / (2 Tc) ds from the root to the latest tip, where k(s) is the number
  of lineages at time s and M the number of mergers (a node with m children counts m - 1). The tip dates are fixed;
  the internal node times are the unknowns. Z(mu) integrates p(k, t | mu) over the internal node times, and

    d log Z / d ln mu = E[sum_b k_b - mu L sum_b tau_b]

  under the posterior of the node times at mu. `--no-prior` drops the Kingman term.

Difference from TreeTime's model
  TreeTime computes the branch likelihoods from the alignment (GTR model, all sites), and its coalescent prior is the
  mean-field approximation with frozen expected lineage counts. This script takes Poisson branch likelihoods from
  realized substitution counts and the coupled Kingman density of the actual node times. Its maximum is therefore the
  maximum of this exact model only, not the target of TreeTime's rate step.

Commands
  profile     maximize log p over the node times for each rate of a grid (a profile likelihood, not Z), then refine
              the joint optimum over rate and node times with Nelder-Mead
  landscape   log p at the true node times and at scaled node depths, and the Hessian at the true times
  mcmc        sample the node times for each rate of a grid in several independent chains; per rate, the pooled
              slope with its Monte Carlo standard error; the zero of the slope over ln mu with its standard error
  validate    run the sampler on cases with known posteriors and compare (see Validation)

Sampler (`mcmc`, `validate`)
  One sweep is a global depth move and one move of every internal node in random order:
  - depth move: all internal node ages before the latest tip are multiplied by f, ln f ~ N(0, s_depth), accepted with
    the Jacobian f^n of the n moved nodes
  - root move: Gaussian random walk with a fixed step s_root
  - other internal nodes: with probability 1/2 a uniform draw between the parent and the earliest child, otherwise a
    Gaussian step with standard deviation 0.2 times that interval
  Every proposal is symmetric or carries its Jacobian, so the chain has the posterior as its stationary law. During
  burn-in, s_depth and s_root adapt toward an acceptance of 0.3; after burn-in they stay fixed, and only post-burn-in
  sweeps enter the estimates. The chains start at the true node times (`--true-dates`) or at times from the branch
  counts at the rate, scaled in depth by the factors of `--start-scales`. Each chain gets its own seed from
  `--seed`, the rate index, and the chain index (numpy SeedSequence); the output records it.

  Monte Carlo standard error: per chain, the sample standard deviation divided by the square root of the effective
  sample size n / tau_int; tau_int is the integrated autocorrelation time from the initial monotone sequence
  estimator of Geyer (1992, https://doi.org/10.1214/ss/1177011137). The variance estimate uses the same estimator on
  (x - mean)^2. Chains of equal length are pooled as the mean of the chain means with
  se = sqrt(sum se_i^2) / n_chains. `chain_p` is the chi-square probability of the spread of the chain means given
  their standard errors; a small value means the chains have not converged to the same law.

  The zero of the slope is interpolated linearly in ln mu between the two grid rates where the pooled slope changes
  sign from positive to negative; its standard error follows from the slope errors by the delta method. The linear
  interpolation adds a bias of second order in the grid spacing; a finer grid around the zero reduces it.

Validation
  `validate` runs 3 chains per case and compares the root age mean and variance and the slope with exact values:
  - cherries with both tips at one date, counts k1, k2 and lambda = mu L: with the Kingman prior the root age is gamma
    with shape k + 1 and rate 2 lambda + 1/Tc, and the slope is k - (k + 1) 2 lambda / (2 lambda + 1/Tc), k = k1 + k2;
    without a prior the rate is 2 lambda and the slope is -1
  - a 3-tip tree ((A,B)X,C)R with tip dates 2025, 2025, 2024.5 and a zero-count branch: moments from a 2-D adaptive
    quadrature of the same density, which checks the interior node moves and the incremental prior update
  A check passes when the sampled value lies within 3 standard errors of the exact value.

  Result of `python test_scripts/poisson_timetree_surface.py validate` (seed 1, 20000 sweeps, 2000 burn-in, 3 chains
  per case): all 21 checks pass (5 cherries, among them k1 = k2 = 0 at two rates and one case without a prior, and 2
  three-tip trees, each with root age mean, root age variance and slope); the largest deviation is 1.40 standard
  errors. A root move whose step depends on the current root age without a Hastings correction fails 17 of the 21
  checks by up to 15.5 standard errors on the same seed.

  Cost: at 400 tips, one sweep takes about 27 ms on one core, and the slope has an autocorrelation time of about 20
  sweeps.

Options
  Shared: --tree <nwk> --dates <tsv with name and date columns> --length <L> [--tc 10 | --no-prior]
          [--true-dates <tsv>] [--rates MIN MAX N] [--seed N]
  mcmc:     [--sweeps 3000] [--burn-in 500] [--start-scales 0.7 1.0 1.4] [--jobs N] [--true-rate MU] [--output <tsv>]
  validate: [--sweeps 20000] [--burn-in 2000] [--seed 1] [--jobs N]

Examples (from the repository root, in the pixi environment)
  python test_scripts/poisson_timetree_surface.py validate
  python test_scripts/poisson_timetree_surface.py mcmc --tree validation/data/const/rep00/L1000_mu0.0003/inputs/realized.nwk \
    --dates validation/data/const/rep00/tip_dates.tsv --true-dates validation/data/const/rep00/node_dates.tsv \
    --length 1000 --tc 10 --rates 1e-4 6e-4 11 --true-rate 3e-4 --seed 1 --jobs 10 --output tmp/mcmc/rep00.tsv
  python test_scripts/poisson_timetree_surface.py profile --tree ... --dates ... --length 1000 --tc 10
"""

import argparse
import io
import logging
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.integrate import quad
from scipy.optimize import LinearConstraint, minimize
from scipy.sparse import coo_matrix
from scipy.special import gammaln
from scipy.stats import chi2

Z_MIN, Z_MAX = -25.0, 6.0
TARGET_ACCEPTANCE = 0.3
PASS_THRESHOLD_SE = 3.0
VALIDATION_TC = 10.0
VALIDATION_LENGTH = 1000
VALIDATION_TIP_DATE = 2025.0

log = logging.getLogger('poisson_timetree_surface')


def main() -> None:
  logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stderr)
  args = parse_args()
  if args.command == 'validate':
    sys.exit(0 if validate(args) else 1)
  tree = TimeTree.from_files(args.tree, args.dates, args.length)
  tc = None if args.no_prior else args.tc
  true_dates = pd.read_csv(args.true_dates, sep='\t', index_col=0)['date'].to_dict() if args.true_dates else None
  log.info(f'{tree.n} nodes, {len(tree.internal)} internal, {int(tree.k.sum())} substitutions on '
           f'{int((tree.k > 0).sum())} of {len(tree.k)} branches; prior: {"none" if tc is None else f"Kingman Tc={tc}"}')
  rates = np.geomspace(args.rates[0], args.rates[1], int(args.rates[2]))
  if args.command == 'mcmc':
    mcmc(tree, tc, true_dates, rates, args)
  elif args.command == 'landscape':
    t_true = tree.true_times(true_dates)
    landscape(tree, t_true, args.tc, rates, args.scales)
    hessian_report(tree, t_true, args.true_rate)
  else:
    profile(tree, tc, true_dates, rates, args.inner)


def parse_args() -> argparse.Namespace:
  p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  commands = p.add_subparsers(dest='command', required=True)
  shared = argparse.ArgumentParser(add_help=False)
  shared.add_argument('--tree', type=Path, required=True)
  shared.add_argument('--dates', type=Path, required=True)
  shared.add_argument('--length', type=int, required=True)
  shared.add_argument('--tc', type=float, default=10.0)
  shared.add_argument('--no-prior', action='store_true')
  shared.add_argument('--true-dates', type=Path)
  shared.add_argument('--rates', type=float, nargs=3, default=[5e-5, 1e-3, 15], metavar=('MIN', 'MAX', 'N'))
  shared.add_argument('--true-rate', type=float, help='true clock rate, to report the maximum relative to it')

  prof = commands.add_parser('profile', parents=[shared], help='profile likelihood over the node times')
  prof.add_argument('--inner', choices=['trust-constr', 'lbfgs', 'nelder-mead'], default='trust-constr')

  land = commands.add_parser('landscape', parents=[shared], help='log p at the true and at scaled node times')
  land.add_argument('--scales', type=float, nargs='+', default=[0.5, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0])

  sampler = argparse.ArgumentParser(add_help=False)
  sampler.add_argument('--seed', type=int, default=1, help='seed of all chains; recorded in the output')
  sampler.add_argument('--jobs', type=int, default=1, help='chains run in parallel')
  sampler.add_argument('--start-scales', type=float, nargs='+', default=[0.7, 1.0, 1.4],
                       help='one chain per factor: start at the base node depths scaled by it')

  chains = commands.add_parser('mcmc', parents=[shared, sampler], help='slope of log Z over ln mu by MCMC')
  chains.add_argument('--sweeps', type=int, default=3000)
  chains.add_argument('--burn-in', type=int, default=500)
  chains.add_argument('--output', type=Path, help='per-rate table as TSV (also printed)')

  check = commands.add_parser('validate', parents=[sampler], help='check the sampler on cases with known posteriors')
  check.add_argument('--sweeps', type=int, default=20000)
  check.add_argument('--burn-in', type=int, default=2000)
  args = p.parse_args()
  if args.command in ('mcmc', 'validate') and args.burn_in >= args.sweeps:
    p.error('--burn-in must be smaller than --sweeps')
  if args.command == 'landscape' and (args.true_dates is None or args.true_rate is None):
    p.error('landscape needs --true-dates and --true-rate')
  return args


class TimeTree:
  def __init__(self, tree: Phylo.BaseTree.Tree, dates: dict[str, float], length: int):
    clades = list(tree.find_clades(order='preorder'))
    index = {id(c): i for i, c in enumerate(clades)}
    self.names = [c.name for c in clades]
    self.n = len(clades)
    self.is_leaf = np.array([c.is_terminal() for c in clades])
    self.parent = np.full(self.n, -1)
    for c in clades:
      for child in c.clades:
        self.parent[index[id(child)]] = index[id(c)]
    self.children = [[index[id(x)] for x in c.clades] for c in clades]
    self.preorder = np.arange(self.n)
    self.postorder = self.preorder[::-1]
    self.internal = np.flatnonzero(~self.is_leaf)
    self.edges = np.flatnonzero(self.parent >= 0)
    self.k = np.array([round((clades[i].branch_length or 0.0) * length) for i in self.edges], dtype=float)
    self.length = length
    self.leaf_times = np.array([dates[self.names[i]] if self.is_leaf[i] else np.nan for i in range(self.n)])
    self.log_k_factorial = gammaln(self.k + 1).sum()
    self.weight = np.where(self.is_leaf, 1.0, 1.0 - np.array([len(c) for c in self.children]))
    self.root = int(np.flatnonzero(self.parent == -1)[0])

  @classmethod
  def from_files(cls, tree_path: Path, dates_path: Path, length: int) -> 'TimeTree':
    dates = pd.read_csv(dates_path, sep='\t', index_col=0)['date'].to_dict()
    return cls(Phylo.read(tree_path, 'newick'), dates, length)

  @classmethod
  def from_newick(cls, newick: str, dates: dict[str, float], length: int) -> 'TimeTree':
    return cls(Phylo.read(io.StringIO(newick), 'newick'), dates, length)

  def true_times(self, true_dates: dict[str, float] | None) -> np.ndarray:
    if true_dates is None:
      raise ValueError('this command needs --true-dates')
    t = self.leaf_times.copy()
    t[self.internal] = [true_dates[self.names[v]] for v in self.internal]
    return t

  def times(self, z):
    t = self.leaf_times.copy()
    argmin = np.full(self.n, -1)
    gap = np.zeros(self.n)
    gap[self.internal] = np.exp(np.clip(z, Z_MIN, Z_MAX))
    for v in self.internal[::-1]:
      kids = self.children[v]
      j = kids[int(np.argmin(t[kids]))]
      argmin[v] = j
      t[v] = t[j] - gap[v]
    parents = self.parent[self.edges]
    tau = (t[self.edges] - t[argmin[parents]]) + gap[parents]
    return t, argmin, tau

  def initial_z(self, mu, t_true=None):
    if t_true is not None:
      gaps = [max(t_true[self.children[v]].min() - t_true[v], 1e-3) for v in self.internal]
      return np.log(gaps)
    k_of = dict(zip(self.edges, self.k))
    return np.log([max(k_of[c] for c in self.children[v]) / (mu * self.length) + 0.05 for v in self.internal])

  def log_lh(self, log_mu, z, tc, with_gradient=False):
    mu = np.exp(log_mu)
    t, argmin, tau = self.times(z)
    rate = mu * self.length
    lh = np.sum(self.k * np.log(rate * tau) - rate * tau) - self.log_k_factorial
    grad_t = np.zeros(self.n)
    if with_gradient:
      d_tau = self.k / tau - rate
      np.add.at(grad_t, self.edges, d_tau)
      np.add.at(grad_t, self.parent[self.edges], -d_tau)
    d_log_mu = np.sum(self.k) - rate * np.sum(tau)
    if tc is not None:
      lh_prior, grad_prior = self.kingman(t, tc, with_gradient)
      lh += lh_prior
      grad_t += grad_prior
    if not with_gradient:
      return lh
    total = grad_t.copy()
    for v in self.preorder:
      if self.is_leaf[v]:
        continue
      j = argmin[v]
      total[j] += total[v]
    order = {v: i for i, v in enumerate(self.internal)}
    inside = (z > Z_MIN) & (z < Z_MAX)
    grad_z = np.array([total[v] * -np.exp(z[order[v]]) for v in self.internal]) * inside
    return lh, d_log_mu, grad_z

  def log_lh_of_times(self, log_mu, t_internal, tc):
    t = self.leaf_times.copy()
    t[self.internal] = t_internal
    tau = t[self.edges] - t[self.parent[self.edges]]
    rate = np.exp(log_mu) * self.length
    with np.errstate(divide='ignore', invalid='ignore'):
      lh = np.sum(np.where(self.k > 0, self.k * np.log(rate * tau), 0.0) - rate * tau) - self.log_k_factorial
      d_tau = np.where(self.k > 0, self.k / tau, 0.0) - rate
    grad_t = np.zeros(self.n)
    np.add.at(grad_t, self.edges, d_tau)
    np.add.at(grad_t, self.parent[self.edges], -d_tau)
    if tc is not None:
      lh_prior, grad_prior = self.kingman(t, tc, True)
      lh += lh_prior
      grad_t += grad_prior
    return lh, grad_t[self.internal]

  def edge_constraints(self):
    position = {v: i for i, v in enumerate(self.internal)}
    rows, cols, vals, lower = [], [], [], []
    for row, c in enumerate(self.edges):
      p = self.parent[c]
      rows.append(row)
      cols.append(position[p])
      vals.append(-1.0)
      if self.is_leaf[c]:
        lower.append(-self.leaf_times[c])
      else:
        rows.append(row)
        cols.append(position[c])
        vals.append(1.0)
        lower.append(0.0)
    matrix = coo_matrix((vals, (rows, cols)), shape=(len(self.edges), len(self.internal))).tocsr()
    return LinearConstraint(matrix, np.array(lower), np.inf)

  def z_of_times(self, t_internal):
    t = self.leaf_times.copy()
    t[self.internal] = t_internal
    return np.log(np.maximum([t[self.children[v]].min() - t[v] for v in self.internal], np.exp(Z_MIN)))

  def earliest_tips(self):
    earliest = self.leaf_times.copy()
    for v in self.internal[::-1]:
      earliest[v] = earliest[self.children[v]].min()
    return earliest

  def scaled_times(self, t, scale):
    anchor = self.earliest_tips()[self.internal] if scale < 1 else np.nanmax(self.leaf_times)
    scaled = t.copy()
    scaled[self.internal] = anchor - scale * (anchor - t[self.internal])
    return scaled

  def hessian(self, log_mu, t):
    tau = t[self.edges] - t[self.parent[self.edges]]
    rate = np.exp(log_mu) * self.length
    position = {v: i + 1 for i, v in enumerate(self.internal)}
    h = np.zeros((len(self.internal) + 1, len(self.internal) + 1))
    for c, k, tb in zip(self.edges, self.k, tau):
      e = np.zeros(len(h))
      if c in position:
        e[position[c]] = 1.0
      e[position[self.parent[c]]] = -1.0
      h -= k / tb**2 * np.outer(e, e)
      h[0] -= rate * e
      h[:, 0] -= rate * e
    h[0, 0] = -rate * tau.sum()
    return h

  def kingman(self, t, tc, with_gradient):
    order = np.argsort(-t, kind='stable')
    ts, ws = t[order], self.weight[order]
    k_after = np.concatenate([[0.0], np.cumsum(ws)[:-1]])
    k_before = k_after + ws
    root = t.min()
    edges = np.append(ts, root)
    durations = edges[:-1] - edges[1:]
    pairs = k_before * (k_before - 1) / 2
    n_mergers = -np.sum(self.weight[~self.is_leaf])
    lh = -np.sum(durations * pairs) / tc - n_mergers * np.log(tc)
    grad = np.zeros(self.n)
    if with_gradient:
      internal = ~self.is_leaf[order]
      grad[order[internal]] = (k_after[internal] * (k_after[internal] - 1) / 2 - pairs[internal]) / tc
    return lh, grad


def profile(tree: TimeTree, tc: float | None, true_dates: dict[str, float] | None, rates: np.ndarray,
            inner: str) -> None:
  rows, z = [], None
  for mu in rates[::-1]:
    z0 = tree.initial_z(mu) if z is None else z
    start = time.time()
    lh, z, nfev = inner_optimum(tree, np.log(mu), z0, tc, inner)
    root, rmse = summarize(tree, z, true_dates)
    rows.append({'mu': mu, 'log_lh': lh, 'T_root': root, 'date_rmse': rmse, 'evals': nfev,
                 'seconds': time.time() - start})
  table = pd.DataFrame(rows[::-1])
  print(table.to_string(index=False, formatters={'mu': '{:.3e}'.format, 'log_lh': '{:.2f}'.format,
                                                 'T_root': '{:.2f}'.format, 'date_rmse': '{:.3f}'.format,
                                                 'seconds': '{:.1f}'.format}))

  best = table.log_lh.idxmax()
  mu0 = table.mu[best]
  _, z_best, _ = inner_optimum(tree, np.log(mu0), tree.initial_z(mu0), tc, inner)
  x0 = np.concatenate([[np.log(mu0)], z_best])
  f = lambda x: -tree.log_lh(x[0], x[1:], tc)
  start = time.time()
  res = minimize(f, x0, method='Nelder-Mead',
                 options={'maxfev': 400 * len(x0), 'xatol': 1e-6, 'fatol': 1e-8, 'adaptive': True})
  root, rmse = summarize(tree, res.x[1:], true_dates)
  print(f'joint optimum (Nelder-Mead from the best grid point, {res.nfev} evaluations, '
        f'{time.time() - start:.0f} s): mu = {np.exp(res.x[0]):.4e}, log L = {-res.fun:.2f}, '
        f'T_root = {root:.2f}, date RMSE = {rmse:.3f}')


def inner_optimum(tree, log_mu, z0, tc, method):
  if method == 'trust-constr':
    t0 = tree.times(z0)[0][tree.internal]
    f = lambda t: tuple(-x for x in tree.log_lh_of_times(log_mu, t, tc))
    res = minimize(f, t0, jac=True, method='trust-constr', constraints=[tree.edge_constraints()],
                   options={'maxiter': 5000, 'gtol': 1e-8, 'xtol': 1e-10})
    return -res.fun, tree.z_of_times(res.x), res.nfev
  if method == 'lbfgs':
    f = lambda z: tuple(-x for x in tree.log_lh(log_mu, z, tc, with_gradient=True)[::2])
    res = minimize(f, np.clip(z0, Z_MIN, Z_MAX), jac=True, method='L-BFGS-B', bounds=[(Z_MIN, Z_MAX)] * len(z0),
                   options={'maxiter': 20000, 'maxfun': 50000, 'gtol': 1e-9, 'ftol': 1e-14})
  else:
    f = lambda z: -tree.log_lh(log_mu, z, tc)
    res = minimize(f, z0, method='Nelder-Mead',
                   options={'maxfev': 400 * len(z0), 'xatol': 1e-6, 'fatol': 1e-8, 'adaptive': True})
  return -res.fun, res.x, res.nfev


def summarize(tree, z, true_dates):
  t, _, _ = tree.times(z)
  root = t[tree.root]
  if true_dates is None:
    return root, np.nan
  truth = np.array([true_dates.get(tree.names[v], np.nan) for v in tree.internal])
  err = t[tree.internal] - truth
  return root, float(np.sqrt(np.nanmean(err**2)))


def landscape(tree, t_true, tc, rates, scales):
  print('log L at the true node times as a function of mu (Poisson only, and with the Kingman prior):')
  grid = []
  for mu in rates:
    poisson = tree.log_lh_of_times(np.log(mu), t_true[tree.internal], None)[0]
    total = tree.log_lh_of_times(np.log(mu), t_true[tree.internal], tc)[0]
    grid.append({'mu': mu, 'poisson': poisson, 'with_prior': total})
  print(pd.DataFrame(grid).to_string(index=False, formatters={'mu': '{:.3e}'.format, 'poisson': '{:.2f}'.format,
                                                              'with_prior': '{:.2f}'.format}))
  for label, prior in [('Poisson only', None), (f'Poisson + Kingman Tc={tc}', tc)]:
    table = {}
    for scale in scales:
      t = tree.scaled_times(t_true, scale)
      tau = t[tree.edges] - t[tree.parent[tree.edges]]
      valid = np.all(tau > 0)
      root = t[tree.root]
      table[f's={scale:g} root={root:.1f}'] = [
        tree.log_lh_of_times(np.log(mu), t[tree.internal], prior)[0] if valid else np.nan for mu in rates
      ]
    df = pd.DataFrame(table, index=[f'{mu:.2e}' for mu in rates]).T
    best = df.stack().idxmax()
    print(f'\n{label}: log L over the scale s of the true node depths (rows) and mu (columns); maximum at {best}')
    print(df.to_string(float_format='{:.1f}'.format))


def hessian_report(tree, t_true, mu):
  h = tree.hessian(np.log(mu), t_true)
  values, vectors = np.linalg.eigh(-h)
  print(f'\nHessian of the Poisson log L in (log mu, {len(tree.internal)} internal times) at the true times, '
        f'mu={mu:.2e}')
  print(f'  eigenvalues of -H: {np.sum(values < 1e-9 * values.max())} zero, smallest nonzero '
        f'{values[values > 1e-9 * values.max()].min():.3g}, largest {values.max():.3g}')
  print(f'  branches with k=0 add no curvature: {int(np.sum(tree.k == 0))} of {len(tree.k)}')
  mu_share = vectors[0] ** 2
  order = np.argsort(-mu_share)[:5]
  print('  eigenvectors with the largest log-mu component:')
  internal_ages = tree.earliest_tips()[tree.internal] - t_true[tree.internal]
  for i in order:
    v = vectors[1:, i]
    corr = np.corrcoef(v, internal_ages)[0, 1] if np.std(v) > 0 else np.nan
    print(f'    eigenvalue {values[i]:.3g}: log-mu share {mu_share[i]:.3f}, '
          f'correlation of the time components with node depth {corr:+.2f}, '
          f'nonzero time components {int(np.sum(np.abs(v) > 1e-3))}')
  return values


def mcmc(tree: TimeTree, tc: float | None, true_dates: dict[str, float] | None, rates: np.ndarray,
         args: argparse.Namespace) -> None:
  settings = ChainSettings(sweeps=args.sweeps, burn_in=args.burn_in)
  jobs = [
    ChainJob(tree, mu, tc, start_times(tree, mu, true_dates, scale), settings, args.seed, (i, j), scale)
    for i, mu in enumerate(rates) for j, scale in enumerate(args.start_scales)
  ]
  log.info(f'{len(jobs)} chains of {settings.sweeps} sweeps on {args.jobs} processes')
  results = run_jobs(jobs, args.jobs)
  rows = [pool_chains(mu, [r for job, r in zip(jobs, results) if job.spawn_key[0] == i]) for i, mu in enumerate(rates)]
  table = pd.DataFrame([row.__dict__ for row in rows]).assign(seed=args.seed)
  if args.output:
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep='\t', index=False, float_format='%.10g')
  print(table.to_string(index=False, float_format='{:.4g}'.format))
  zero = slope_zero(np.log(rates), table.slope.to_numpy(), table.slope_se.to_numpy())
  if zero is None:
    log.warning('the pooled slope does not change sign from + to - on the rate grid: no maximum inside the grid')
    return
  ratio = f', {zero.mu / args.true_rate:.4f} times the true rate {args.true_rate:g}' if args.true_rate else ''
  print(f'\nmaximum of Z: ln mu = {zero.ln_mu:.5f} +- {zero.ln_mu_se:.5f} (Monte Carlo SE), '
        f'mu = {zero.mu:.5e}{ratio}')
  if zero.crossings > 1:
    log.warning(f'the slope changes sign {zero.crossings} times on the grid; reported the first crossing')


def start_times(tree: TimeTree, mu: float, true_dates: dict[str, float] | None, scale: float) -> np.ndarray:
  base = tree.true_times(true_dates) if true_dates is not None else tree.times(tree.initial_z(mu))[0]
  return tree.scaled_times(base, scale)


@dataclass(frozen=True)
class ChainSettings:
  sweeps: int
  burn_in: int
  thin: int = 5


@dataclass(frozen=True)
class ChainJob:
  tree: TimeTree
  mu: float
  tc: float | None
  t_start: np.ndarray
  settings: ChainSettings
  seed: int
  spawn_key: tuple[int, ...]
  start_scale: float


def run_jobs(jobs: list[ChainJob], n_jobs: int) -> list['ChainResult']:
  if n_jobs <= 1:
    return [run_chain(job) for job in jobs]
  with ProcessPoolExecutor(n_jobs) as pool:
    return list(pool.map(run_chain, jobs))


def run_chain(job: ChainJob) -> 'ChainResult':
  started = time.time()
  rng = np.random.default_rng(np.random.SeedSequence(job.seed, spawn_key=job.spawn_key))
  sampler = Sampler(job.tree, job.mu, job.tc, job.t_start)
  slope, root, samples = [], [], []
  for sweep in range(job.settings.sweeps):
    adapt = sweep < job.settings.burn_in
    sampler.sweep(rng, sweep + 1 if adapt else None)
    if adapt:
      continue
    slope.append(sampler.slope_term())
    root.append(sampler.t[job.tree.root])
    if sweep % job.settings.thin == 0:
      samples.append(sampler.t[job.tree.internal].copy())
  slope, root = np.array(slope), np.array(root)
  root_var = (root - root.mean()) ** 2
  return ChainResult(
    mu=job.mu, start_scale=job.start_scale, spawn_key=job.spawn_key,
    slope=slope.mean(), slope_se=mcse(slope), slope_ess=effective_sample_size(slope),
    root_mean=root.mean(), root_mean_se=mcse(root), root_var=root_var.mean(), root_var_se=mcse(root_var),
    acceptance=sampler.node_accepted / max(sampler.node_proposed, 1),
    root_acceptance=sampler.root_accepted / max(sampler.root_proposed, 1),
    depth_acceptance=sampler.depth_accepted / max(sampler.depth_proposed, 1),
    root_step=sampler.root_step, depth_step=sampler.depth_step,
    covariance_term=covariance_term(job.tree, np.array(samples), job.tc),
    seconds=time.time() - started,
  )


@dataclass(frozen=True)
class ChainResult:
  mu: float
  start_scale: float
  spawn_key: tuple[int, ...]
  slope: float
  slope_se: float
  slope_ess: float
  root_mean: float
  root_mean_se: float
  root_var: float
  root_var_se: float
  acceptance: float
  root_acceptance: float
  depth_acceptance: float
  root_step: float
  depth_step: float
  covariance_term: float
  seconds: float


class Sampler:
  """Metropolis-Hastings chain over the internal node times at a fixed rate."""

  def __init__(self, tree: TimeTree, mu: float, tc: float | None, t_start: np.ndarray):
    tau = t_start[tree.edges] - t_start[tree.parent[tree.edges]]
    if np.any(tau <= 0):
      raise ValueError('the start times put a node at or after one of its children')
    self.tree = tree
    self.tc = tc
    self.rate = mu * tree.length
    self.t = t_start.copy()
    self.k_of = np.zeros(tree.n)
    self.k_of[tree.edges] = tree.k
    self.latest = np.nanmax(tree.leaf_times)
    self.order = tree.internal.copy()
    self.root_step = max(float(np.std(t_start[tree.internal])), 1.0)
    self.depth_step = 0.02
    self.node_accepted = self.node_proposed = 0
    self.root_accepted = self.root_proposed = 0
    self.depth_accepted = self.depth_proposed = 0

  def sweep(self, rng: np.random.Generator, adapt_iteration: int | None) -> None:
    accepted = self._depth_move(rng)
    if adapt_iteration is not None:
      self.depth_step = adapted(self.depth_step, accepted, adapt_iteration)
    rng.shuffle(self.order)
    for v in self.order:
      if v == self.tree.root:
        accepted = self._root_move(rng)
        if adapt_iteration is not None:
          self.root_step = adapted(self.root_step, accepted, adapt_iteration)
      else:
        self._node_move(rng, v)

  def slope_term(self) -> float:
    tau = self.t[self.tree.edges] - self.t[self.tree.parent[self.tree.edges]]
    return self.tree.k.sum() - self.rate * tau.sum()

  def _depth_move(self, rng: np.random.Generator) -> bool:
    self.depth_proposed += 1
    log_f = rng.normal(0.0, self.depth_step)
    proposal = self.t.copy()
    internal = self.tree.internal
    proposal[internal] = self.latest - np.exp(log_f) * (self.latest - self.t[internal])
    delta = self._log_target(proposal) - self._log_target(self.t) + len(internal) * log_f
    if np.log(rng.random()) < delta:
      self.t = proposal
      self.depth_accepted += 1
      return True
    return False

  def _root_move(self, rng: np.random.Generator) -> bool:
    self.root_proposed += 1
    v = self.tree.root
    new = self.t[v] + rng.normal(0.0, self.root_step)
    return self._accept_node(rng, v, new, upper=min(self.t[c] for c in self.tree.children[v]))

  def _node_move(self, rng: np.random.Generator, v: int) -> None:
    self.node_proposed += 1
    upper = min(self.t[c] for c in self.tree.children[v])
    lower = self.t[self.tree.parent[v]]
    if rng.random() < 0.5:
      new = rng.uniform(lower, upper)
    else:
      new = self.t[v] + rng.normal(0.0, 0.2 * (upper - lower))
    if not new > lower:
      return
    self.node_accepted += self._accept_node(rng, v, new, upper)

  def _accept_node(self, rng: np.random.Generator, v: int, new: float, upper: float) -> bool:
    if not new < upper:
      return False
    old = self.t[v]
    delta = self._local_log_lh(v, new) - self._local_log_lh(v, old)
    delta += self._kingman_change(v, old, new) if self.tc is not None else 0.0
    if np.log(rng.random()) < delta:
      self.t[v] = new
      if v == self.tree.root:
        self.root_accepted += 1
      return True
    return False

  def _local_log_lh(self, v: int, x: float) -> float:
    tree = self.tree
    branches = ([(self.k_of[v], x - self.t[tree.parent[v]])] if v != tree.root else []) + [
      (self.k_of[c], self.t[c] - x) for c in tree.children[v]
    ]
    return sum(poisson_log_term(k, self.rate * tau) for k, tau in branches)

  def _kingman_change(self, v: int, old: float, new: float) -> float:
    """Change of log p_Kingman when node v moves from old to new; the lineage count changes by its weight on [a, b]."""
    tree = self.tree
    a, b = min(old, new), max(old, new)
    others = np.arange(tree.n) != v
    k_integral = np.sum(tree.weight[others] * np.clip(self.t[others] - a, 0.0, b - a))
    w = tree.weight[v]
    integral = w * k_integral + 0.5 * w * (w - 1) * (b - a)
    return -integral / self.tc if new > old else integral / self.tc

  def _log_target(self, times: np.ndarray) -> float:
    tree = self.tree
    tau = times[tree.edges] - times[tree.parent[tree.edges]]
    if np.any(tau <= 0):
      return -np.inf
    poisson = np.sum(np.where(tree.k > 0, tree.k * np.log(self.rate * tau), 0.0) - self.rate * tau)
    return poisson + (tree.kingman(times, self.tc, False)[0] if self.tc is not None else 0.0)


def poisson_log_term(k: float, expected: float) -> float:
  return k * np.log(expected) - expected if k > 0 else -expected


def adapted(step: float, accepted: bool, iteration: int) -> float:
  return step * np.exp((float(accepted) - TARGET_ACCEPTANCE) / np.sqrt(iteration))


def covariance_term(tree: TimeTree, samples: np.ndarray, tc: float | None) -> float:
  """Part of E[int k(k-1)/(2 Tc)] from the covariance of node times, which a product of node marginals leaves out."""
  if tc is None or len(samples) < 2:
    return np.nan
  latest = np.nanmax(tree.leaf_times)
  grid = np.linspace(samples.min() - 1.0, latest, 1000)
  leaf_part = np.sum(tree.leaf_times[tree.is_leaf][None, :] > grid[:, None], axis=1)
  w_int = tree.weight[tree.internal]
  above = samples[:, :, None] > grid[None, None, :]
  k = leaf_part[None, :] + np.einsum('i,sig->sg', w_int, above)
  p_surv = above.mean(axis=0)
  own_var = np.einsum('i,ig->g', w_int**2, p_surv * (1 - p_surv))
  return 0.5 * np.trapezoid(k.var(axis=0) - own_var, grid) / tc


def effective_sample_size(x: np.ndarray) -> float:
  """n / tau_int, with tau_int from the initial monotone sequence estimator (Geyer 1992)."""
  n = len(x)
  centered = x - x.mean()
  variance = centered @ centered / n
  if variance == 0:
    return float(n)
  size = 1 << (2 * n - 1).bit_length()
  spectrum = np.fft.rfft(centered, size)
  rho = np.fft.irfft(spectrum * np.conj(spectrum), size)[:n] / (n * variance)
  pairs = rho[: (n // 2) * 2].reshape(-1, 2).sum(axis=1)
  positive = np.flatnonzero(pairs <= 0)
  pairs = pairs[: positive[0] if len(positive) else len(pairs)]
  tau = -1.0 + 2.0 * np.minimum.accumulate(pairs).sum()
  return n / max(tau, 1.0 / np.log10(max(n, 10)))


def mcse(x: np.ndarray) -> float:
  return float(np.std(x) / np.sqrt(effective_sample_size(x)))


def pool_chains(mu: float, chains: list[ChainResult]) -> 'PooledRate':
  means = np.array([c.slope for c in chains])
  ses = np.array([c.slope_se for c in chains])
  pooled = means.mean()
  dispersion = np.sum(((means - pooled) / ses) ** 2)
  return PooledRate(
    mu=mu,
    slope=pooled,
    slope_se=np.sqrt(np.sum(ses**2)) / len(chains),
    chain_slopes=','.join(f'{m:.4g}' for m in means),
    chain_p=float(chi2.sf(dispersion, len(chains) - 1)) if len(chains) > 1 else np.nan,
    slope_ess=sum(c.slope_ess for c in chains),
    root_mean=np.mean([c.root_mean for c in chains]),
    root_sd=np.sqrt(np.mean([c.root_var for c in chains])),
    acceptance=np.mean([c.acceptance for c in chains]),
    root_acceptance=np.mean([c.root_acceptance for c in chains]),
    depth_acceptance=np.mean([c.depth_acceptance for c in chains]),
    covariance_term=np.mean([c.covariance_term for c in chains]),
    seconds=sum(c.seconds for c in chains),
  )


@dataclass(frozen=True)
class PooledRate:
  mu: float
  slope: float
  slope_se: float
  chain_slopes: str
  chain_p: float
  slope_ess: float
  root_mean: float
  root_sd: float
  acceptance: float
  root_acceptance: float
  depth_acceptance: float
  covariance_term: float
  seconds: float


def slope_zero(ln_mu: np.ndarray, slope: np.ndarray, se: np.ndarray) -> 'SlopeZero | None':
  crossings = np.flatnonzero((slope[:-1] > 0) & (slope[1:] <= 0))
  if len(crossings) == 0:
    return None
  i = crossings[0]
  x1, x2, s1, s2 = ln_mu[i], ln_mu[i + 1], slope[i], slope[i + 1]
  x0 = x1 + (x2 - x1) * s1 / (s1 - s2)
  variance = (x2 - x1) ** 2 * (s2**2 * se[i] ** 2 + s1**2 * se[i + 1] ** 2) / (s1 - s2) ** 4
  return SlopeZero(ln_mu=x0, ln_mu_se=float(np.sqrt(variance)), mu=float(np.exp(x0)), crossings=len(crossings))


@dataclass(frozen=True)
class SlopeZero:
  ln_mu: float
  ln_mu_se: float
  mu: float
  crossings: int


def validate(args: argparse.Namespace) -> bool:
  cases = validation_cases()
  settings = ChainSettings(sweeps=args.sweeps, burn_in=args.burn_in)
  jobs = [
    ChainJob(case.tree, case.mu, case.tc, case.tree.scaled_times(case.t_start, scale), settings, args.seed, (i, j),
             scale)
    for i, case in enumerate(cases) for j, scale in enumerate(args.start_scales)
  ]
  results = run_jobs(jobs, args.jobs)
  rows = []
  for i, case in enumerate(cases):
    chains = [r for job, r in zip(jobs, results) if job.spawn_key[0] == i]
    sampled = {
      'root age mean': (VALIDATION_TIP_DATE - np.mean([c.root_mean for c in chains]),
                        pooled_se([c.root_mean_se for c in chains])),
      'root age variance': (np.mean([c.root_var for c in chains]), pooled_se([c.root_var_se for c in chains])),
      'slope': (np.mean([c.slope for c in chains]), pooled_se([c.slope_se for c in chains])),
    }
    for quantity, exact in case.exact.items():
      value, se = sampled[quantity]
      deviation = (value - exact) / se
      rows.append({'case': case.name, 'quantity': quantity, 'exact': exact, 'sampled': value, 'se': se,
                   'deviation_se': deviation, 'pass': abs(deviation) < PASS_THRESHOLD_SE})
  table = pd.DataFrame(rows)
  print(f'seed {args.seed}, {args.sweeps} sweeps, {args.burn_in} burn-in, {len(args.start_scales)} chains per case')
  print(table.to_string(index=False, float_format='{:.6g}'.format))
  failed = int((~table['pass']).sum())
  if failed:
    log.error(f'{failed} of {len(table)} checks deviate by {PASS_THRESHOLD_SE:g} standard errors or more')
  return failed == 0


def pooled_se(ses: list[float]) -> float:
  return float(np.sqrt(np.sum(np.square(ses))) / len(ses))


@dataclass(frozen=True)
class ValidationCase:
  name: str
  tree: TimeTree
  mu: float
  tc: float | None
  t_start: np.ndarray
  exact: dict[str, float]


def validation_cases() -> list[ValidationCase]:
  cherries = [(0, 0, 3e-4, VALIDATION_TC), (0, 0, 1e-3, VALIDATION_TC), (2, 1, 3e-4, VALIDATION_TC),
              (12, 9, 1e-3, VALIDATION_TC), (3, 4, 3e-4, None)]
  three_tip = [(0, 3, 2, 1, 3e-4, VALIDATION_TC), (4, 2, 0, 5, 1e-3, VALIDATION_TC)]
  return [cherry_case(*c) for c in cherries] + [three_tip_case(*c) for c in three_tip]


def cherry_case(k1: int, k2: int, mu: float, tc: float | None) -> ValidationCase:
  length = VALIDATION_LENGTH
  tree = TimeTree.from_newick(f'(A:{k1 / length},B:{k2 / length})R;', {'A': VALIDATION_TIP_DATE,
                                                                       'B': VALIDATION_TIP_DATE}, length)
  lam = mu * length
  k = k1 + k2
  rate = 2 * lam + (1 / tc if tc is not None else 0.0)
  t_start = tree.leaf_times.copy()
  t_start[tree.root] = VALIDATION_TIP_DATE - (k + 1) / rate
  exact = {'root age mean': (k + 1) / rate, 'root age variance': (k + 1) / rate**2,
           'slope': k - (k + 1) * 2 * lam / rate}
  prior = f'Tc={tc:g}' if tc is not None else 'no prior'
  return ValidationCase(f'cherry k=({k1},{k2}) mu={mu:g} {prior}', tree, mu, tc, t_start, exact)


def three_tip_case(k_a: int, k_b: int, k_x: int, k_c: int, mu: float, tc: float) -> ValidationCase:
  length = VALIDATION_LENGTH
  date_c = VALIDATION_TIP_DATE - 0.5
  newick = f'((A:{k_a / length},B:{k_b / length})X:{k_x / length},C:{k_c / length})R;'
  dates = {'A': VALIDATION_TIP_DATE, 'B': VALIDATION_TIP_DATE, 'C': date_c}
  tree = TimeTree.from_newick(newick, dates, length)
  lam = mu * length
  counts = {'A': k_a, 'B': k_b, 'X': k_x, 'C': k_c}

  def log_density(x: float, r: float) -> float:
    taus = {'A': VALIDATION_TIP_DATE - x, 'B': VALIDATION_TIP_DATE - x, 'X': x - r, 'C': date_c - r}
    poisson = sum(poisson_log_term(counts[b], lam * tau) for b, tau in taus.items())
    return poisson + kingman_log_density_reference([VALIDATION_TIP_DATE, VALIDATION_TIP_DATE, date_c], [x, r], tc)

  def slope(x: float, r: float) -> float:
    total_tau = 2 * (VALIDATION_TIP_DATE - x) + (x - r) + (date_c - r)
    return sum(counts.values()) - lam * total_tau

  exact = posterior_moments_2d(log_density, slope, date_c)
  t_start = tree.leaf_times.copy()
  t_start[tree.names.index('X')] = exact['x mean']
  t_start[tree.root] = min(exact['x mean'], date_c) - 1.0
  moments = {'root age mean': exact['root age mean'], 'root age variance': exact['root age variance'],
             'slope': exact['slope']}
  return ValidationCase(f'3 tips k=(A{k_a},B{k_b},X{k_x},C{k_c}) mu={mu:g} Tc={tc:g}', tree, mu, tc, t_start,
                        moments)


def kingman_log_density_reference(tip_dates: list[float], merge_dates: list[float], tc: float) -> float:
  """Kingman log density of binary mergers, summed interval by interval from the latest tip back to the root."""
  events = sorted([(d, +1) for d in tip_dates] + [(d, -1) for d in merge_dates], key=lambda e: -e[0])
  lineages, integral = 0, 0.0
  for (date, change), (next_date, _) in zip(events, events[1:]):
    lineages += change
    integral += lineages * (lineages - 1) / 2 * (date - next_date)
  return -len(merge_dates) * np.log(tc) - integral / tc


def posterior_moments_2d(log_density, slope, date_c: float) -> dict[str, float]:
  """Moments of the posterior over (t_X, t_R) by nested adaptive quadrature.

  The outer integral over t_X is split at the date of C, and the inner one runs over the root age below the earliest
  child, so each piece has a smooth integrand; the unbounded pieces use the infinite-range rule of QUADPACK.
  The integrand is scaled by its peak, so the absolute error bound of 1e-9 is small against the integrals.
  """
  peak = max(log_density(x, min(x, date_c) - u) for x in np.linspace(date_c - 20, VALIDATION_TIP_DATE - 1e-6, 100)
             for u in np.geomspace(1e-3, 50, 100))

  def inner(x: float, g) -> float:
    upper = min(x, date_c)
    return checked_quad(lambda u: g(x, upper - u) * np.exp(log_density(x, upper - u) - peak), 0, np.inf)

  def integral(g) -> float:
    recent = checked_quad(lambda x: inner(x, g), date_c, VALIDATION_TIP_DATE)
    old = checked_quad(lambda v: inner(date_c - v, g), 0, np.inf)
    return recent + old

  z = integral(lambda x, r: 1.0)
  age_mean = integral(lambda x, r: VALIDATION_TIP_DATE - r) / z
  age_sq = integral(lambda x, r: (VALIDATION_TIP_DATE - r - age_mean) ** 2) / z
  return {
    'root age mean': age_mean,
    'root age variance': age_sq,
    'slope': integral(slope) / z,
    'x mean': integral(lambda x, r: x) / z,
  }


def checked_quad(f, a: float, b: float) -> float:
  value, error = quad(f, a, b, epsabs=1e-12, epsrel=1e-10, limit=200)
  if error > 1e-9 + 1e-8 * abs(value):
    raise ArithmeticError(f'quadrature on [{a}, {b}] reached an error of {error:.2g} at the value {value:.6g}')
  return value


if __name__ == '__main__':
  main()
