"""Shape of TreeTime's log marginal likelihood of node times over the clock rate, and the three-probe rate step.

Purpose
  The marginal rate step of timetree (`fn step_marginal_clock_rate()`) fits a parabola through log Z at the rate and
  at the rate divided and multiplied by 1.25. This script measures, on profiles of log Z over the rate, how quadratic
  log Z is in four coordinates u of the rate (rate, sqrt, cube root, log), and where one parabola step, and repeated
  unguarded full steps, land in each coordinate.

Input
  A TSV with columns L, mu, seed, x, log_z, where x = ln(rate / true rate), as printed by the ignored test
  `research_rate_profile_shape` of `packages/treetime/src/test_utils/__tests__/test_rate_profile.rs` (branch
  `research/rate-step`): log Z under the coalescent prior frozen at the true node times, with the branch likelihoods
  at the true rate (`FrozenRateState`).

Output
  Per (L, mu), medians over seeds: the maximum x_max, the curvature in x, the standardized third derivative per
  coordinate (0 for a parabola), RMS residuals of quadratic fits over windows around the maximum, and the error and
  convergence of the steps. Writes <input>_shape.tsv and <input>_steps.tsv.

Usage
  pixi run python test_scripts/timetree_rate_step_shape.py profile.tsv
"""

import sys

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

H = np.log(1.25)
COORDS = {'rate': 1.0, 'sqrt': 0.5, 'cbrt': 1 / 3, 'log': 0.0}  # u = exp(p x), p = 0: u = x


def to_u(x, p):
    return x if p == 0 else np.exp(p * x)


def to_x(u, p):
    if p == 0:
        return u
    return np.log(u) / p if u > 0 else np.nan


def vertex(f, x0, p):
    """Vertex of the parabola in u through the probes at x0 - H, x0, x0 + H; nan without a maximum."""
    xs = np.array([x0 - H, x0, x0 + H])
    u, y = to_u(xs, p), f(xs)
    a, b, _ = np.polyfit(u, y, 2)
    return to_x(-b / (2 * a), p) if a < 0 else np.nan


def iterate(f, x0, p, lo, hi, n=20, tol=0.01):
    """Unguarded full parabola steps; steps to |dx| < tol, or nan if it leaves the profile or finds no maximum."""
    x = x0
    for i in range(1, n + 1):
        xn = vertex(f, x, p)
        if not np.isfinite(xn) or not lo + H <= xn <= hi - H:
            return np.nan, xn
        if abs(xn - x) < tol:
            return i, xn
        x = xn
    return np.nan, x


def main():
    d = pd.read_csv(sys.argv[1], sep='\t', comment='#')
    d = d[pd.to_numeric(d['L'], errors='coerce').notna()].astype(float)
    shape, steps = [], []
    for (L, mu, seed), g in d.groupby(['L', 'mu', 'seed']):
        g = g.sort_values('x')
        g = g[np.isfinite(g.log_z)]
        f = CubicSpline(g.x.values, g.log_z.values)
        lo, hi = g.x.min(), g.x.max()
        opt = minimize_scalar(lambda x: -f(x), bounds=(lo, hi), method='bounded')
        xs = opt.x
        interior = lo + 0.3 < xs < hi - 0.3
        row = {'L': L, 'mu': mu, 'seed': seed, 'x_max': xs, 'curv_x': float(f(xs, 2))}
        for name, p in COORDS.items():
            # standardized third derivative in u at the maximum, from derivatives in x
            if p == 0:
                d2, d3 = f(xs, 2), f(xs, 3)
            else:  # u = e^{px}: d/du = (1/(p u)) d/dx
                u = np.exp(p * xs)
                f1, f2, f3 = f(xs, 1), f(xs, 2), f(xs, 3)
                d2 = (f2 - p * f1) / (p * u) ** 2
                d3 = (f3 - 3 * p * f2 + 2 * p**2 * f1) / (p * u) ** 3
            row[f'skew_{name}'] = float(d3 / (-d2) ** 1.5) if d2 < 0 else np.nan
            for w in (0.3, 0.7, 1.2):
                m = (g.x >= xs - w) & (g.x <= xs + w)
                if m.sum() >= 5 and xs - w >= lo and xs + w <= hi:
                    u = to_u(g.x[m].values, p)
                    res = g.log_z[m].values - np.polyval(np.polyfit(u, g.log_z[m].values, 2), u)
                    row[f'rms_{name}_w{w}'] = float(np.sqrt(np.mean(res**2)))
        shape.append(row)
        if not interior:
            continue
        for dx in (-1.4, -0.7, 0.7, 1.4):
            x0 = xs + dx
            if not lo + H <= x0 <= hi - H:
                continue
            for name, p in COORDS.items():
                x1 = vertex(f, x0, p)
                n, xf = iterate(f, x0, p, lo, hi)
                steps.append({'L': L, 'mu': mu, 'seed': seed, 'start': dx, 'coord': name,
                              'err_1step': x1 - xs if np.isfinite(x1) else np.nan, 'n_steps': n,
                              'err_final': xf - xs if np.isfinite(n) else np.nan})  # fmt: skip
    shape, steps = pd.DataFrame(shape), pd.DataFrame(steps)
    pd.set_option('display.width', 250)
    print('## maximum of log Z (x = ln rate/true), curvature in x, standardized 3rd derivative per coordinate')
    cols = ['x_max', 'curv_x'] + [f'skew_{c}' for c in COORDS]
    print(shape.groupby(['L', 'mu'])[cols].median().round(3).to_string())
    print('\n## RMS residual (nats) of a quadratic fit in u over |x - x_max| <= w, median over seeds')
    for w in (0.3, 0.7, 1.2):
        c = [f'rms_{n}_w{w}' for n in COORDS]
        print(f'w = {w}')
        print(shape.groupby(['L', 'mu'])[c].median().rename(columns=dict(zip(c, COORDS))).round(3).to_string())
    if len(steps):
        print('\n## one parabola step: |error| in x after the step, median over seeds; failures (no maximum) count')
        t = steps.assign(abs_err=steps.err_1step.abs(), fail=steps.err_1step.isna())
        print(t.pivot_table(index=['L', 'mu', 'start'], columns='coord', values='abs_err', aggfunc='median').round(3).to_string())
        print('\nno maximum after one step (count)')
        print(t.pivot_table(index=['L', 'mu', 'start'], columns='coord', values='fail', aggfunc='sum').to_string())
        print('\n## unguarded full steps: median steps to converge (nan: diverged or no maximum), and diverged count')
        print(steps.pivot_table(index=['L', 'mu', 'start'], columns='coord', values='n_steps', aggfunc='median').to_string())
        print(steps.assign(div=steps.n_steps.isna()).pivot_table(index=['L', 'mu', 'start'], columns='coord', values='div', aggfunc='sum').to_string())
    shape.to_csv(sys.argv[1].replace('.tsv', '_shape.tsv'), sep='\t', index=False)
    steps.to_csv(sys.argv[1].replace('.tsv', '_steps.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
