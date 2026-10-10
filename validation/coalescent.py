"""Kingman coalescent with heterochronous sampling and time-varying coalescent timescale Tc(tau).

tau is time before the most recent sample (years). A pair of lineages merges at rate 1/Tc(tau).
Each schedule provides the cumulative hazard Lambda(tau) = int_0^tau 1/Tc and its inverse.
"""

import numpy as np
from Bio.Phylo.BaseTree import Clade, Tree


class Constant:
    def __init__(self, tc=10.0):
        self.tc = tc

    def cum(self, tau):
        return tau / self.tc

    def inv(self, lam):
        return lam * self.tc


class Exponential:
    """Tc(tau) = tc0 * exp(-tau / scale)."""

    def __init__(self, tc0=10.0, scale=5.0):
        self.tc0, self.scale = tc0, scale

    def cum(self, tau):
        return self.scale / self.tc0 * np.expm1(tau / self.scale)

    def inv(self, lam):
        return self.scale * np.log1p(lam * self.tc0 / self.scale)


class Piecewise:
    """Tc = values[i] for tau in [breaks[i-1], breaks[i]); the last value extends to infinity."""

    def __init__(self, breaks=(5.0, 6.0), values=(10.0, 1.0, 3.0)):
        self.breaks = np.array([0.0, *breaks])
        self.values = np.array(values)
        self.lam = np.concatenate([[0.0], np.cumsum(np.diff(self.breaks) / self.values[:-1])])

    def cum(self, tau):
        i = np.searchsorted(self.breaks, tau, side='right') - 1
        return self.lam[i] + (tau - self.breaks[i]) / self.values[i]

    def inv(self, lam):
        i = np.searchsorted(self.lam, lam, side='right') - 1
        return self.breaks[i] + (lam - self.lam[i]) * self.values[i]


SCHEDULES = {'const': Constant, 'exp': Exponential, 'piecewise': Piecewise}


def simulate(schedule, rng, n_per_batch=100, n_batches=40, interval=0.1):
    """Return (tree, tau) where tau maps every clade to its time before the most recent sample."""
    sample_times = [i * interval for i in range(n_batches)]
    tau, active, n_tips, n_nodes = {}, [], 0, 0
    t, next_batch = 0.0, 0

    while True:
        t_sample = sample_times[next_batch] if next_batch < n_batches else np.inf
        k = len(active)
        t_coal = np.inf
        if k > 1:
            t_coal = schedule.inv(schedule.cum(t) + rng.exponential() / (k * (k - 1) / 2))
        if t_sample <= t_coal:
            for _ in range(n_per_batch):
                n_tips += 1
                c = Clade(name=f'tip_{n_tips:04d}')
                tau[c] = t_sample
                active.append(c)
            t, next_batch = t_sample, next_batch + 1
            continue
        t = t_coal
        i, j = sorted(rng.choice(k, size=2, replace=False), reverse=True)
        children = [active.pop(i), active.pop(j)]
        n_nodes += 1
        parent = Clade(name=f'node_{n_nodes:04d}', clades=children)
        tau[parent] = t
        for c in children:
            c.branch_length = t - tau[c]
        active.append(parent)
        if len(active) == 1 and next_batch == n_batches:
            break

    root = active[0]
    root.name, root.branch_length = 'node_root', 0.0
    return Tree(root=root, rooted=True), tau
