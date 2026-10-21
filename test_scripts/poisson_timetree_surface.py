"""Joint likelihood of the clock rate and the internal node times of a fixed tree, with Poisson branch likelihoods.

log L(mu, t) = sum_b [k_b ln(mu L tau_b) - mu L tau_b - ln k_b!] + log p_Kingman(t | Tc)

k_b are the substitutions on branch b (branch length times L), tau_b = t_child - t_parent, and the Kingman density
of heterochronous samples is sum_mergers (-ln Tc) - int k(k-1)/(2 Tc) dt from the root to the latest tip.
Internal node times are parametrized as t_v = min(children) - exp(z_v), so every branch stays positive; z is bounded
to [Z_MIN, Z_MAX], and each branch length is computed from the gaps so that it does not cancel against calendar years.

For each rate on a grid, the node times are optimized (the profile likelihood, a maximum and not an integral over
node times), and the joint optimum over rate and node times is refined with Nelder-Mead.

Usage: python3 poisson_timetree_surface.py --tree realized.nwk --dates tip_dates.tsv --length 1000 --tc 10 \
         [--true-dates node_dates.tsv] [--no-prior] [--rates 5e-5 1e-3 15] [--inner lbfgs|nelder-mead]
"""

import argparse
import time

import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.optimize import LinearConstraint, minimize
from scipy.sparse import coo_matrix
from scipy.special import gammaln

Z_MIN, Z_MAX = -25.0, 6.0


class TimeTree:
  def __init__(self, tree_path, dates_path, length):
    tree = Phylo.read(tree_path, 'newick')
    dates = pd.read_csv(dates_path, sep='\t', index_col=0)['date']
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
    h[0, 0] = -rate * tau.sum()
    for c, k, tb in zip(self.edges, self.k, tau):
      e = np.zeros(len(h))
      if c in position:
        e[position[c]] = 1.0
      e[position[self.parent[c]]] = -1.0
      h -= k / tb**2 * np.outer(e, e)
      h[0] -= rate * e
      h[:, 0] -= rate * e
      h[0, 0] += 0.0
    h[0, 0] = -rate * tau.sum()
    return h

  def kingman(self, t, tc, with_gradient):
    weight = np.where(self.is_leaf, 1.0, 1.0 - np.array([len(c) for c in self.children]))
    order = np.argsort(-t, kind='stable')
    ts, ws = t[order], weight[order]
    k_after = np.concatenate([[0.0], np.cumsum(ws)[:-1]])
    k_before = k_after + ws
    root = t.min()
    edges = np.append(ts, root)
    durations = edges[:-1] - edges[1:]
    pairs = k_before * (k_before - 1) / 2
    n_mergers = -np.sum(weight[~self.is_leaf])
    lh = -np.sum(durations * pairs) / tc - n_mergers * np.log(tc)
    grad = np.zeros(self.n)
    if with_gradient:
      internal = ~self.is_leaf[order]
      grad[order[internal]] = (k_after[internal] * (k_after[internal] - 1) / 2 - pairs[internal]) / tc
    return lh, grad


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
      root = t[tree.parent == -1][0]
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
  print(f'\nHessian of the Poisson log L in (log mu, 399 internal times) at the true times, mu={mu:.2e}')
  print(f'  eigenvalues of -H: {np.sum(values < 1e-9 * values.max())} zero, smallest nonzero '
        f'{values[values > 1e-9 * values.max()].min():.3g}, largest {values.max():.3g}')
  tau = t_true[tree.edges] - t_true[tree.parent[tree.edges]]
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


def mcmc_slope(tree, t_start, mu, tc, sweeps, burn_in, seed):
  rng = np.random.default_rng(seed)
  t = t_start.copy()
  rate = mu * tree.length
  weight = np.where(tree.is_leaf, 1.0, 1.0 - np.array([len(c) for c in tree.children]))
  k_of = np.zeros(tree.n)
  k_of[tree.edges] = tree.k
  is_root = tree.parent == -1

  def local_log_lh(v, x):
    total = 0.0
    if not is_root[v]:
      tau = x - t[tree.parent[v]]
      total += k_of[v] * np.log(rate * tau) - rate * tau if k_of[v] > 0 else -rate * tau
    for c in tree.children[v]:
      tau = t[c] - x
      total += k_of[c] * np.log(rate * tau) - rate * tau if k_of[c] > 0 else -rate * tau
    return total

  def kingman_change(v, old, new):
    a, b = min(old, new), max(old, new)
    others = np.arange(tree.n) != v
    k_integral = np.sum(weight[others] * np.clip(t[others] - a, 0.0, b - a))
    integral = weight[v] * k_integral + 0.5 * weight[v] * (weight[v] - 1) * (b - a)
    return -integral / tc if new > old else integral / tc

  def full_log_target(times):
    tau = times[tree.edges] - times[tree.parent[tree.edges]]
    if np.any(tau <= 0):
      return -np.inf
    poisson = np.sum(np.where(tree.k > 0, tree.k * np.log(rate * np.maximum(tau, 1e-300)), 0.0) - rate * tau)
    return poisson + (tree.kingman(times, tc, False)[0] if tc is not None else 0.0)

  latest = np.nanmax(tree.leaf_times)
  samples = []
  record = []
  accepted = 0
  scale_accepted = 0
  order = tree.internal.copy()
  for sweep in range(sweeps):
    log_f = rng.normal(0.0, 0.02)
    proposal = t.copy()
    proposal[tree.internal] = latest - np.exp(log_f) * (latest - t[tree.internal])
    delta = full_log_target(proposal) - full_log_target(t) + len(tree.internal) * log_f
    if np.log(rng.random()) < delta:
      t = proposal
      scale_accepted += 1
    rng.shuffle(order)
    for v in order:
      upper = min(t[c] for c in tree.children[v])
      lower = -np.inf if is_root[v] else t[tree.parent[v]]
      old = t[v]
      if np.isfinite(lower) and rng.random() < 0.5:
        new = rng.uniform(lower, upper)
      else:
        width = (upper - lower) if np.isfinite(lower) else max(upper - old, 0.5)
        new = old + rng.normal(0.0, 0.2 * width)
      if not lower < new < upper:
        continue
      delta = local_log_lh(v, new) - local_log_lh(v, old)
      delta += kingman_change(v, old, new) if tc is not None else 0.0
      if np.log(rng.random()) < delta:
        t[v] = new
        accepted += 1
    if sweep >= burn_in:
      tau = t[tree.edges] - t[tree.parent[tree.edges]]
      record.append((tree.k.sum() - rate * tau.sum(), t[is_root][0]))
      if sweep % 5 == 0:
        samples.append(t[tree.internal].copy())
  samples = np.array(samples)
  grid = np.linspace(samples.min() - 1.0, latest, 1000)
  leaf_part = np.sum(tree.leaf_times[tree.is_leaf][None, :] > grid[:, None], axis=1)
  w_int = weight[tree.internal]
  above = samples[:, :, None] > grid[None, None, :]
  k = leaf_part[None, :] + np.einsum('i,sig->sg', w_int, above)
  p_surv = above.mean(axis=0)
  var_k = k.var(axis=0)
  own_var = np.einsum('i,ig->g', w_int**2, p_surv * (1 - p_surv))
  covariance_term = 0.5 * np.trapz(var_k - own_var, grid) / tc if tc is not None else np.nan
  record = np.array(record)
  batches = np.array_split(record[:, 0], 20)
  means = np.array([b.mean() for b in batches])
  return {
    'mu': mu,
    'slope': record[:, 0].mean(),
    'slope_se': means.std(ddof=1) / np.sqrt(len(means)),
    'root_mean': record[:, 1].mean(),
    'root_sd': record[:, 1].std(),
    'acceptance': accepted / (sweeps * len(tree.internal)),
    'scale_acceptance': scale_accepted / sweeps,
    'covariance_term': covariance_term,
  }


def summarize(tree, z, true_dates):
  t, _, _ = tree.times(z)
  root = t[tree.parent == -1][0]
  if true_dates is None:
    return root, np.nan
  truth = np.array([true_dates.get(tree.names[v], np.nan) for v in tree.internal])
  err = t[tree.internal] - truth
  return root, float(np.sqrt(np.nanmean(err**2)))


def main():
  p = argparse.ArgumentParser()
  p.add_argument('--tree', required=True)
  p.add_argument('--dates', required=True)
  p.add_argument('--length', type=int, required=True)
  p.add_argument('--tc', type=float, default=10.0)
  p.add_argument('--no-prior', action='store_true')
  p.add_argument('--true-dates')
  p.add_argument('--rates', type=float, nargs=3, default=[5e-5, 1e-3, 15], metavar=('MIN', 'MAX', 'N'))
  p.add_argument('--inner', choices=['trust-constr', 'lbfgs', 'nelder-mead'], default='trust-constr')
  p.add_argument('--landscape', action='store_true', help='evaluate at the true node times instead of optimizing')
  p.add_argument('--scales', type=float, nargs='+', default=[0.5, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0])
  p.add_argument('--true-rate', type=float, default=3e-4)
  p.add_argument('--mcmc', action='store_true', help='d log Z / d ln mu = E[sum k - mu L sum tau] by MCMC over node times')
  p.add_argument('--sweeps', type=int, default=3000)
  p.add_argument('--start-scale', type=float, default=1.0, help='start the chain at the true depths scaled by this')
  p.add_argument('--burn-in', type=int, default=500)
  args = p.parse_args()

  tree = TimeTree(args.tree, args.dates, args.length)
  tc = None if args.no_prior else args.tc
  true_dates = pd.read_csv(args.true_dates, sep='\t', index_col=0)['date'].to_dict() if args.true_dates else None
  print(f'{tree.n} nodes, {len(tree.internal)} internal, {int(tree.k.sum())} substitutions on '
        f'{int((tree.k > 0).sum())} of {len(tree.k)} branches; prior: {"none" if tc is None else f"Kingman Tc={tc}"}')

  rates = np.geomspace(args.rates[0], args.rates[1], int(args.rates[2]))
  if args.mcmc:
    t_true = tree.leaf_times.copy()
    t_true[tree.internal] = [true_dates[tree.names[v]] for v in tree.internal]
    t_start = tree.scaled_times(t_true, args.start_scale)
    rows = [mcmc_slope(tree, t_start, mu, tc, args.sweeps, args.burn_in, seed=i) for i, mu in enumerate(rates)]
    print(pd.DataFrame(rows).to_string(index=False, formatters={'mu': '{:.3e}'.format, 'slope': '{:.2f}'.format,
                                                                'slope_se': '{:.2f}'.format, 'root_mean': '{:.2f}'.format,
                                                                'root_sd': '{:.2f}'.format, 'acceptance': '{:.3f}'.format,
                                                                'scale_acceptance': '{:.3f}'.format,
                                                                'covariance_term': '{:.2f}'.format}))
    return
  if args.landscape:
    t_true = tree.leaf_times.copy()
    t_true[tree.internal] = [true_dates[tree.names[v]] for v in tree.internal]
    landscape(tree, t_true, args.tc, rates, args.scales)
    hessian_report(tree, t_true, args.true_rate)
    return
  rows, z = [], None
  for mu in rates[::-1]:
    z0 = tree.initial_z(mu) if z is None else z
    start = time.time()
    lh, z, nfev = inner_optimum(tree, np.log(mu), z0, tc, args.inner)
    root, rmse = summarize(tree, z, true_dates)
    rows.append({'mu': mu, 'log_lh': lh, 'T_root': root, 'date_rmse': rmse, 'evals': nfev,
                 'seconds': time.time() - start})
  profile = pd.DataFrame(rows[::-1])
  print(profile.to_string(index=False, formatters={'mu': '{:.3e}'.format, 'log_lh': '{:.2f}'.format,
                                                   'T_root': '{:.2f}'.format, 'date_rmse': '{:.3f}'.format,
                                                   'seconds': '{:.1f}'.format}))

  best = profile.log_lh.idxmax()
  mu0 = profile.mu[best]
  _, z_best, _ = inner_optimum(tree, np.log(mu0), tree.initial_z(mu0), tc, args.inner)
  x0 = np.concatenate([[np.log(mu0)], z_best])
  f = lambda x: -tree.log_lh(x[0], x[1:], tc)
  start = time.time()
  res = minimize(f, x0, method='Nelder-Mead',
                 options={'maxfev': 400 * len(x0), 'xatol': 1e-6, 'fatol': 1e-8, 'adaptive': True})
  root, rmse = summarize(tree, res.x[1:], true_dates)
  print(f'joint optimum (Nelder-Mead from the best grid point, {res.nfev} evaluations, '
        f'{time.time() - start:.0f} s): mu = {np.exp(res.x[0]):.4e}, log L = {-res.fun:.2f}, '
        f'T_root = {root:.2f}, date RMSE = {rmse:.3f}')


if __name__ == '__main__':
  main()
