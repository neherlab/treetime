"""Summary tables and figures for the validation results, optionally compared against a baseline run.

Usage: python analyze.py [--results results] [--baseline other/results] [--out results/report]
Reads results/{ancestral,timetree}.tsv written by evaluate.py; writes summary.md, compare.md and *.png to --out.
"""

import argparse
import itertools
import warnings
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from coalescent import PRESENT, SCHEDULES  # noqa: E402

HERE = Path(__file__).parent
CORE = ['-'.join(v) for v in itertools.product(['keep', 'reroot'], ['none', 'opt', 'skyline'], ['est', 'fixed'])]
COALS = ['const', 'exp', 'piecewise']

# Categorical slots 1-3 of the reference palette (validated light; aqua is below 3:1, so every series also has its
# own marker and a legend). Colors follow the entity: a prior keeps its color in every figure.
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
PRIOR_STYLE = {'none': (BLUE, 'o'), 'opt': (ORANGE, 's'), 'skyline': (AQUA, '^')}
ANC_STYLE = {'marginal': (BLUE, 'o'), 'marginal-gtr3': (ORANGE, 's'), 'parsimony': (AQUA, '^')}
L_STYLE = {1000: (BLUE, 'o'), 10000: (ORANGE, 's')}
INK, INK2, GRID, SURFACE = '#0b0b0b', '#52514e', '#e6e5e0', '#fcfcfb'
# Ordinal blue ramp (steps 250/400/550/700, validated --ordinal): substitution rate, light = low.
MU_COLOR = {0.0003: '#86b6ef', 0.001: '#3987e5', 0.003: '#1c5cab', 0.01: '#0d366b'}
L_LINESTYLE = {1000: '--', 10000: '-'}
TRUTH = '#c3c2b7'

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE, 'font.size': 9,
    'text.color': INK, 'axes.labelcolor': INK2, 'xtick.color': INK2, 'ytick.color': INK2, 'axes.edgecolor': GRID,
    'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
    'lines.linewidth': 1.5, 'lines.markersize': 6, 'legend.frameon': False, 'axes.titlesize': 9.5,
})  # fmt: skip


def load(results):
    tables = {}
    for suite in ['ancestral', 'timetree', 'optimize']:
        path = results / f'{suite}.tsv'
        if path.exists():
            d = pd.read_csv(path, sep='\t')
            if suite == 'timetree':
                core = d.variant.isin(CORE)
                d[['root', 'prior', 'clock']] = d.variant.where(core, 'extra-extra-extra').str.split('-', expand=True)
            tables[suite] = d
    return tables


def md(df, floatfmt='.3g'):
    """Markdown table with row keys as columns and column keys as `name=value` labels."""
    df = df.copy()
    names = df.columns.names
    if any(n is not None for n in names):
        keys = df.columns if isinstance(df.columns, pd.MultiIndex) else [(c,) for c in df.columns]
        df.columns = [' '.join(f'{n}={v}' for n, v in zip(names, k)) for k in keys]
    return df.reset_index().to_markdown(index=False, floatfmt=floatfmt) + '\n'


def q10_50_90(x):
    return f'{x.median():.2f} [{x.quantile(0.1):.2f}, {x.quantile(0.9):.2f}]' if len(x) else ''


def summary(t):
    out = ['# Validation summary\n']
    for suite, d in t.items():
        out.append(f'**{suite}**: {len(d)} runs, {(d.returncode != 0).sum()} failed\n')
    for suite, d in t.items():
        f = d[d.returncode != 0]
        if len(f):
            out += [f'\n## Failures: {suite}\n', md(f.groupby(['error', 'coal', 'L']).size().unstack('L', fill_value=0))]

    if 'ancestral' in t:
        a = t['ancestral'][t['ancestral'].returncode == 0]
        out += ['\n## Ancestral\n', 'Per-site error of internal sequences (median over replicates)\n',
                md(a.pivot_table(index='variant', columns=['L', 'mu'], values='site_error', aggfunc='median')),
                '\nInferred GTR (median; mu is in units of the input branch lengths, truth 1)\n',
                md(a[a.variant != 'parsimony'].groupby(['variant', 'L', 'mu'])[
                    ['pi_max_abs_error', 'W_log_rms_error', 'mu_inferred']].median())]  # fmt: skip

    if 'timetree' in t:
        d = t['timetree'][t['timetree'].returncode == 0]
        c, keep = d[d.variant.isin(CORE)], d[d.variant.isin(CORE) & (d.root == 'keep')]
        out += ['\n## Timetree\n',
                f'Nodes dated before their parent: {int(d.n_negative_branches.sum())}; fraction of inferred nodes '
                f'matching a true clade: min {keep.frac_nodes_matched.min():.4f} (keep root), '
                f'{c[c.root == "reroot"].frac_nodes_matched.min():.4f} (reroot)\n',
                '\nInternal-node date RMSE in years (keep root, median over replicates)\n',
                md(keep.pivot_table(index=['clock', 'prior', 'coal'], columns=['L', 'mu'], values='date_rmse', aggfunc='median')),
                '\nInternal-node date bias in years (keep root, median)\n',
                md(keep.pivot_table(index=['clock', 'prior'], columns=['L', 'mu'], values='date_bias', aggfunc='median')),
                '\nEstimated clock rate, relative error: median [10%, 90%] (keep root, opt)\n',
                md(keep[(keep.clock == 'est') & (keep.prior == 'opt')].pivot_table(
                    index=['coal', 'L'], columns='mu', values='clock_rate_rel_error', aggfunc=q10_50_90)),
                '\nTrue root recovered by `--reroot` (fraction of runs)\n',
                md(c[c.root == 'reroot'].pivot_table(index='coal', columns=['L', 'mu'], values='root_recovered', aggfunc='mean')),
                '\nCoalescent Tc: log bias, mean |log error|, band coverage (keep root, fixed clock; means)\n',
                md(keep[(keep.prior != 'none') & (keep.clock == 'fixed')].groupby(['prior', 'coal', 'L'])[
                    ['tc_log_bias', 'tc_log_mae', 'tc_coverage']].mean())]  # fmt: skip
        x = d[~d.variant.isin(CORE)]
        if len(x):
            out += ['\nExtra variants: date RMSE (median), 90% CI coverage and width (mean)\n',
                    md(x.groupby(['variant', 'mu'])[['date_rmse', 'ci90_coverage', 'ci90_mean_width']].agg(
                        {'date_rmse': 'median', 'ci90_coverage': 'mean', 'ci90_mean_width': 'mean'}))]  # fmt: skip
        out += ['\nRun time in seconds and peak memory in MB (median, max)\n',
                md(t['timetree'].groupby(['L', 'mu'])[['seconds', 'max_rss_mb']].agg(['median', 'max']))]
        bad = c[(c.L * c.mu >= 3) & (c.date_rmse > 1)]
        if len(bad):
            out += ['\nOutliers: date RMSE > 1 y although L*mu >= 3\n',
                    md(bad.groupby(['coal', 'rep', 'L', 'mu', 'clock']).agg(
                        runs=('variant', 'size'), date_rmse=('date_rmse', 'median'),
                        clock_rate_rel_error=('clock_rate_rel_error', 'median')))]  # fmt: skip
    if 'optimize' in t:
        o = t['optimize'][t['optimize'].returncode == 0]
        fixed = o[o.variant == 'fixed-topo']
        out += ['\n## Branch length optimization\n',
                'Per-branch RMSE in substitutions against the realized branch lengths (median over replicates; '
                'the two root branches excluded). `input` is the perturbed starting tree\n',
                md(o.pivot_table(index=['variant', 'coal'], columns=['L', 'mu'], values='rmse_subs', aggfunc='median')),
                md(fixed.pivot_table(index='coal', columns=['L', 'mu'], values='input_rmse_subs', aggfunc='median')
                   .rename(index=lambda c: f'input, {c}')),
                '\nTotal length ratio, optimized over realized and over expected (mu t); fixed topology, median\n',
                md(fixed.groupby(['L', 'mu'])[['length_ratio', 'length_ratio_expected']].median()),
                '\nZero branches kept at zero, non-zero branches kept non-zero, and squared deviation from mu t relative '
                'to Poisson noise (about 1 when optimal); fixed topology, mean\n',
                md(fixed.groupby(['L', 'mu'])[['zero_recall', 'nonzero_recall', 'dispersion_vs_expected']].mean()),
                '\nSum of the two root branches, relative error against realized and expected: median [10%, 90%]; '
                'fixed topology\n',
                md(fixed.groupby(['L', 'mu']).agg(
                    realized=('root_sum_rel_error', q10_50_90), expected=('root_sum_rel_error_expected', q10_50_90))),
                '\nFraction of true branches present after optimization with default settings (mean)\n',
                md(o[o.variant == 'default'].pivot_table(index='coal', columns=['L', 'mu'],
                                                         values='frac_branches_matched', aggfunc='mean'))]  # fmt: skip
    return ''.join(out)


COMPARE = {
    'ancestral': ['site_error', 'W_log_rms_error', 'mu_inferred', 'seconds'],
    'timetree': ['date_rmse', 'date_bias', 'clock_rate_rel_error', 'tc_log_mae', 'seconds'],
    'optimize': ['rmse_subs', 'length_ratio', 'root_sum_rel_error', 'seconds'],
}
KEY = ['variant', 'coal', 'rep', 'L', 'mu']


def compare(t, base):
    out = ['# Comparison against baseline\n', 'Differences are current minus baseline, matched per run.\n']
    for suite, d in t.items():
        if suite not in base:
            continue
        m = d.merge(base[suite], on=KEY, suffixes=('', '_base'))
        now_ok, was_ok = m.returncode == 0, m.returncode_base == 0
        out.append(f'\n## {suite}: {len(m)} matched runs; newly failing {(~now_ok & was_ok).sum()}, '
                   f'newly passing {(now_ok & ~was_ok).sum()}\n')  # fmt: skip
        if (~now_ok & was_ok).any():
            out += ['\nNewly failing\n', md(m[~now_ok & was_ok].groupby(['error', 'variant']).size().to_frame('runs'))]
        both = m[now_ok & was_ok]
        diff = pd.DataFrame({c: both[c] - both[f'{c}_base'] for c in COMPARE[suite] if c != 'seconds'})
        diff['seconds_ratio'] = both.seconds / both.seconds_base
        diff[['variant', 'L', 'mu']] = both[['variant', 'L', 'mu']]
        agg = {c: 'median' for c in diff if c not in ('variant', 'L', 'mu')}
        out += ['\nMedian change by variant\n', md(diff.groupby('variant').agg(agg))]
    return ''.join(out)


def legend_handles(styles, label_fmt='{}'):
    return [plt.Line2D([], [], color=c, marker=m, label=label_fmt.format(k)) for k, (c, m) in styles.items()]


def fig_date_rmse(d, path):
    keep = d[(d.returncode == 0) & d.variant.isin(CORE) & (d.root == 'keep')]
    fig, axes = plt.subplots(2, 3, figsize=(10, 6), sharex=True, sharey=True)
    for (i, L), (j, coal) in itertools.product(enumerate([1000, 10000]), enumerate(COALS)):
        ax = axes[i, j]
        for (prior, (color, marker)), (clock, ls) in itertools.product(PRIOR_STYLE.items(), [('fixed', '-'), ('est', '--')]):
            s = keep[(keep.L == L) & (keep.coal == coal) & (keep.prior == prior) & (keep.clock == clock)]
            med = s.groupby('mu').date_rmse.median()
            ax.plot(med.index, med.values, ls=ls, color=color, marker=marker, mfc=SURFACE if clock == 'est' else color)
        ax.set(xscale='log', yscale='log', title=f'{coal}, L={L}')
        if j == 0:
            ax.set_ylabel('internal-node date RMSE (years)')
        if i == 1:
            ax.set_xlabel('substitution rate mu (per site and year)')
    handles = legend_handles(PRIOR_STYLE, 'prior: {}') + [
        plt.Line2D([], [], color=INK2, ls='-', label='clock fixed at true mu'),
        plt.Line2D([], [], color=INK2, ls='--', marker='o', mfc=SURFACE, label='clock estimated')]  # fmt: skip
    fig.legend(handles=handles, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle('Dating accuracy (true root kept; median over replicates)', y=1.05, x=0.06, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def fig_clock_rate(d, path):
    s = d[(d.returncode == 0) & (d.variant == 'keep-opt-est')]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4), sharey=True)
    mus = sorted(s.mu.unique())
    for ax, coal in zip(axes, COALS):
        for k, (L, (color, marker)) in enumerate(L_STYLE.items()):
            for x, mu in enumerate(mus):
                r = 1 + s[(s.coal == coal) & (s.L == L) & (s.mu == mu)].clock_rate_rel_error
                xs = x + (k - 0.5) * 0.3 + np.linspace(-0.06, 0.06, len(r))
                ax.scatter(xs, r, color=color, marker=marker, s=22, edgecolor=SURFACE, linewidth=0.8, zorder=3)
                ax.plot([x + (k - 0.5) * 0.3 - 0.1, x + (k - 0.5) * 0.3 + 0.1], [r.median()] * 2, color=INK, lw=1.5)
        ax.axhline(1, color=INK2, lw=1, zorder=1)
        ax.set(yscale='log', title=coal, xticks=range(len(mus)), xticklabels=[f'{m:g}' for m in mus])
        ax.set_xlabel('substitution rate mu')
    axes[0].set_ylabel('estimated / true clock rate')
    fig.legend(handles=legend_handles(L_STYLE, 'L={}') + [plt.Line2D([], [], color=INK, label='median')],
               loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.02))  # fmt: skip
    fig.suptitle('Estimated clock rate (keep root, opt prior; failed runs absent)', y=1.1, x=0.06, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def fig_tc(d, path):
    s = d[(d.returncode == 0) & d.variant.isin(CORE) & (d.root == 'keep') & (d.clock == 'fixed') & (d.prior != 'none')]
    fig, axes = plt.subplots(2, 3, figsize=(10, 5.5), sharex=True, sharey=True)
    for (i, L), (j, coal) in itertools.product(enumerate([1000, 10000]), enumerate(COALS)):
        ax = axes[i, j]
        for k, prior in enumerate(['opt', 'skyline']):
            color, marker = PRIOR_STYLE[prior]
            g = s[(s.L == L) & (s.coal == coal) & (s.prior == prior)].groupby('mu').tc_log_bias
            x = g.median().index * (1.08 if k else 1 / 1.08)
            ax.errorbar(x, g.median(), yerr=[g.median() - g.quantile(0.25), g.quantile(0.75) - g.median()],
                        color=color, marker=marker, capsize=0, elinewidth=1)  # fmt: skip
        ax.axhline(0, color=INK2, lw=1)
        ax.set(xscale='log', title=f'{coal}, L={L}')
        if j == 0:
            ax.set_ylabel('log(Tc inferred / Tc true)')
        if i == 1:
            ax.set_xlabel('substitution rate mu')
    fig.legend(handles=legend_handles({p: PRIOR_STYLE[p] for p in ['opt', 'skyline']}, 'prior: {}'),
               loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.0))  # fmt: skip
    fig.suptitle('Coalescent time scale (keep root, clock fixed; median and interquartile range)', y=1.05, x=0.06, ha='left')
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def fig_ancestral(a, path):
    a = a[a.returncode == 0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for (variant, (color, marker)), (L, ls) in itertools.product(ANC_STYLE.items(), [(1000, '--'), (10000, '-')]):
        s = a[(a.variant == variant) & (a.L == L)].groupby('mu')
        mfc = SURFACE if L == 1000 else color
        axes[0].plot(s.site_error.median().index, s.site_error.median(), ls=ls, color=color, marker=marker, mfc=mfc)
        if variant != 'parsimony':
            axes[1].plot(s.mu_inferred.median().index, s.mu_inferred.median(), ls=ls, color=color, marker=marker, mfc=mfc)
    axes[0].set(xscale='log', yscale='log', xlabel='substitution rate mu', ylabel='per-site error, internal nodes',
                title='Reconstruction error (median)')  # fmt: skip
    axes[1].axhline(1, color=INK2, lw=1)
    axes[1].set(xscale='log', yscale='log', xlabel='substitution rate mu', ylabel='inferred GTR mu (truth 1)',
                title='Inferred rate scale (median)')  # fmt: skip
    handles = legend_handles(ANC_STYLE) + [
        plt.Line2D([], [], color=INK2, ls='-', label='L=10000'),
        plt.Line2D([], [], color=INK2, ls='--', marker='o', mfc=SURFACE, label='L=1000')]  # fmt: skip
    fig.legend(handles=handles, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 1.04))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def fig_skylines(results, coal, path, variant='keep-skyline-fixed'):
    """Median inferred skyline per (L, mu) against the true Tc(t); colour follows mu, line style follows L."""
    runs = sorted((results / 'timetree' / variant / coal).glob('rep*/L*/coalescent.tsv'))
    if not runs:
        return
    sky = {r: pd.read_csv(r, sep='\t') for r in runs}
    t = np.linspace(min(s['segment.start'].min() for s in sky.values()), PRESENT, 500)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot(t, SCHEDULES[coal]().tc(PRESENT - t), color=TRUTH, lw=7, solid_capstyle='butt', zorder=1)
    for L, ls in L_LINESTYLE.items():
        for mu, color in MU_COLOR.items():
            curves = []
            for r, s in sky.items():
                if r.parent.name != f'L{L}_mu{mu:g}':
                    continue
                i = np.searchsorted(s['segment.start'].values, t, side='right') - 1
                inside = (i >= 0) & (t <= s['segment.end'].values[-1])
                curves.append(np.where(inside, s['T_c.value'].values[np.clip(i, 0, None)], np.nan))
            if not curves:
                continue
            curves = np.array(curves)
            enough = np.sum(~np.isnan(curves), axis=0) >= len(curves) / 2  # at least half the replicates reach t
            with warnings.catch_warnings():  # dates no replicate reaches are all-NaN columns
                warnings.simplefilter('ignore', RuntimeWarning)
                med = np.where(enough, np.nanmedian(curves, axis=0), np.nan)
            ax.plot(t, med, color=color, ls=ls, lw=1.5, zorder=2)
    ax.set(yscale='log', xlabel='date', ylabel='Tc (years)')
    handles = ([plt.Line2D([], [], color=TRUTH, lw=7, label='true Tc')]
               + [plt.Line2D([], [], color=c, lw=1.5, label=f'mu={mu:g}') for mu, c in MU_COLOR.items()]
               + [plt.Line2D([], [], color=INK2, lw=1.5, ls=ls, label=f'L={L}') for L, ls in L_LINESTYLE.items()])  # fmt: skip
    ax.legend(handles=handles, loc='center left', bbox_to_anchor=(1.01, 0.5))
    ax.set_title(f'Skyline against truth: {coal} (true root kept, clock fixed; median over replicates)', loc='left')
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def slope_guide(ax, x, y, label='slope 1/2'):
    """Dotted line of slope 1/2 on log-log axes through the point (x, y)."""
    xs = np.geomspace(x / 30, x, 20)
    ax.plot(xs, y * np.sqrt(xs / x), color=INK2, ls=':', lw=1.5, zorder=1)
    ax.annotate(label, (xs[8], y * np.sqrt(xs[8] / x)), textcoords='offset points', xytext=(-10, 8), color=INK2)


def fig_optimize(o, path):
    """Branch length optimization with the true topology, against mu L (proportional to the substitution count)."""
    o = o[(o.returncode == 0) & (o.variant == 'fixed-topo')].assign(muL=lambda d: d.mu * d.L)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    coal_style = {'const': (BLUE, 'o'), 'exp': (ORANGE, 's'), 'piecewise': (AQUA, '^')}
    panels = [('rmse_subs_expected', 'poisson_rmse_subs', 'RMSE of optimized - mu t (substitutions)',
               'Against expected mu t (grey: Poisson expectation)'),
              ('rmse_subs', 'input_rmse_subs', 'RMSE of optimized - realized (substitutions)',
               'Against realized differences (grey: perturbed input)')]  # fmt: skip
    for ax, (value, reference, ylabel, title) in zip(axes, panels):
        for (coal, (color, marker)), (L, ls) in itertools.product(coal_style.items(), L_LINESTYLE.items()):
            s = o[(o.coal == coal) & (o.L == L)].groupby('muL')
            mfc = SURFACE if L == 1000 else color
            ax.plot(s[value].median().index, s[value].median(), color=color, marker=marker, ls=ls, mfc=mfc)
            ax.plot(s[reference].median().index, s[reference].median(), color=GRID, ls=ls, lw=1.5, zorder=1)
        ax.set(xscale='log', yscale='log', xlabel='mu L', ylabel=ylabel, title=title)
    top = o.groupby('muL').rmse_subs_expected.median()
    slope_guide(axes[0], top.index[-1], top.iloc[-1] * 2)
    top = o.groupby('muL').rmse_subs.median()
    slope_guide(axes[1], top.index[-1], top.iloc[-1] * 2)

    for k, (L, (color, marker)) in enumerate(L_STYLE.items()):
        for muL, g in o[o.L == L].groupby('muL'):
            r = g.root_sum_rel_error_expected.dropna()
            xs = muL * np.exp((k - 0.5) * 0.25 + np.linspace(-0.07, 0.07, len(r)))
            axes[2].scatter(xs, r, color=color, marker=marker, s=14, edgecolor=SURFACE, linewidth=0.6, zorder=3)
            xm = muL * np.exp((k - 0.5) * 0.25)
            axes[2].plot([xm / 1.1, xm * 1.1], [r.median()] * 2, color=INK, lw=1.5)
    axes[2].axhline(0, color=INK2, lw=1)
    axes[2].set(xscale='log', xlabel='mu L', ylabel='relative error of the sum (one dot per replicate)',
                title='Sum of the two root branches vs mu t')  # fmt: skip
    axes[2].legend(handles=[plt.Line2D([], [], color=c, marker=m, ls='', label=f'L={L}') for L, (c, m) in L_STYLE.items()]
                   + [plt.Line2D([], [], color=INK, label='median')], loc='upper right')  # fmt: skip
    handles = legend_handles(coal_style) + [
        plt.Line2D([], [], color=INK2, ls=ls, label=f'L={L}') for L, ls in L_LINESTYLE.items()]  # fmt: skip
    fig.legend(handles=handles, loc='upper center', ncol=5, bbox_to_anchor=(0.36, 1.04))
    fig.suptitle('Branch length optimization (true topology kept). RMSE over the branches of one tree, median over '
                 'replicates; mu L is proportional to the number of substitutions', y=1.1, x=0.04, ha='left')  # fmt: skip
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def fig_runtime(d, path):
    s = d[d.returncode == 0]
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for L, (color, marker) in L_STYLE.items():
        g = s[s.L == L].groupby('mu').seconds
        ax.errorbar(g.median().index * L, g.median(), yerr=[g.median() - g.quantile(0.1), g.quantile(0.9) - g.median()],
                    color=color, marker=marker, capsize=0, elinewidth=1, label=f'L={L}')  # fmt: skip
    ax.set(xscale='log', yscale='log', xlabel='L x mu (expected substitutions per year, whole sequence)',
           ylabel='timetree run time (s)', title='Run time: median and 10-90% range over all variants')  # fmt: skip
    ax.legend(loc='upper left')
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--results', type=Path, default=HERE / 'results')
    p.add_argument('--baseline', type=Path, help='results directory of an earlier run to compare against')
    p.add_argument('--out', type=Path, help='default: <results>/report')
    args = p.parse_args()
    out = args.out or args.results / 'report'
    out.mkdir(parents=True, exist_ok=True)

    t = load(args.results)
    (out / 'summary.md').write_text(summary(t))
    if 'timetree' in t:
        fig_date_rmse(t['timetree'], out / 'timetree_date_rmse.png')
        fig_clock_rate(t['timetree'], out / 'timetree_clock_rate.png')
        fig_tc(t['timetree'], out / 'timetree_tc.png')
        fig_runtime(t['timetree'], out / 'timetree_runtime.png')
        for coal in COALS:
            fig_skylines(args.results, coal, out / f'timetree_skyline_{coal}.png')
    if 'ancestral' in t:
        fig_ancestral(t['ancestral'], out / 'ancestral.png')
    if 'optimize' in t:
        fig_optimize(t['optimize'], out / 'optimize.png')
    if args.baseline:
        (out / 'compare.md').write_text(compare(t, load(args.baseline)))
    print(f'wrote {out}')


if __name__ == '__main__':
    main()
