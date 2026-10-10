"""Compute metrics for every finished run and write results/<suite>.tsv; print a summary table.

Usage: python evaluate.py --suite timetree [-j 16]
Per-run metrics are cached in metrics.json next to run.json.
"""

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

from coalescent import PRESENT, SCHEDULES
from inputs import read_fasta_xz
from run import available_cpus

HERE = Path(__file__).parent
AGE_BINS = [0, 1, 3, 10, np.inf]  # true node age in years before the most recent sample


def age_bin_labels():
    return [f'{a:g}-{b:g}y' for a, b in zip(AGE_BINS[:-1], AGE_BINS[1:])]


def tip_sets(tree):
    sets = {}
    for c in tree.find_clades(order='postorder'):
        sets[c] = frozenset([c.name]) if c.is_terminal() else frozenset().union(*(sets[x] for x in c))
    return sets


def true_dates(rep_dir):
    tables = [pd.read_csv(rep_dir / f, sep='\t', index_col=0) for f in ['node_dates.tsv', 'tip_dates.tsv']]
    return pd.concat(tables)


def binned(err, tau, prefix):
    out = {}
    for label, (a, b) in zip(age_bin_labels(), zip(AGE_BINS[:-1], AGE_BINS[1:])):
        e = err[(tau >= a) & (tau < b)]
        out[f'{prefix}_{label}'] = float(np.sqrt(np.mean(e**2))) if len(e) else np.nan
    return out


def eval_ancestral(out, set_dir):
    rep_dir = set_dir.parent
    truth, inferred = read_fasta_xz(set_dir / 'internal.fasta.xz'), read_fasta_xz(out / 'anc.fasta.xz')
    tau = true_dates(rep_dir)['tau']
    names = sorted(truth)
    err = pd.Series({n: np.mean(truth[n] != inferred[n]) for n in names})
    m = {'n_nodes': len(names), 'n_missing': len(set(names) - set(inferred)), 'site_error': float(err.mean())}
    for label, (a, b) in zip(age_bin_labels(), zip(AGE_BINS[:-1], AGE_BINS[1:])):
        e = err[(tau[names] >= a) & (tau[names] < b)]
        m[f'site_error_{label}'] = float(e.mean()) if len(e) else np.nan

    if not (out / 'gtr.json').exists():
        return m
    g_true, g_inf = (json.loads(p.read_text()) for p in [set_dir / 'gtr.json', out / 'gtr.json'])
    W_true, W_inf = np.array(g_true['W']), np.array(g_inf['W'])
    iu = np.triu_indices(4, 1)
    m['pi_max_abs_error'] = float(np.max(np.abs(np.array(g_inf['pi']) - np.array(g_true['pi']))))
    m['W_log_rms_error'] = float(np.sqrt(np.mean(np.log(W_inf[iu] / W_true[iu]) ** 2)))
    m['mu_inferred'] = g_inf['mu']  # branch lengths are substitutions per site, so the truth is 1
    return m


def eval_timetree(out, set_dir, coal):
    rep_dir = set_dir.parent
    mu = json.loads((set_dir / 'gtr.json').read_text())['mu']
    dates = true_dates(rep_dir)
    true_tree = Phylo.read(rep_dir / 'tree.nwk', 'newick')
    true_by_set = {s: c.name for c, s in tip_sets(true_tree).items() if not c.is_terminal()}

    tree = Phylo.read(out / 'timetree.nwk', 'newick')
    nodes = json.loads((out / 'node_data.json').read_text())['nodes']
    sets = tip_sets(tree)
    m = {'clock_rate_rel_error': json.loads((out / 'clock.json').read_text())['clock_rate'] / mu - 1}

    numdate = {c: nodes[c.name]['numdate'] for c in tree.find_clades()}
    parent = {c: p for p in tree.get_nonterminals() for c in p}
    m['n_negative_branches'] = sum(numdate[c] < numdate[p] - 1e-9 for c, p in parent.items())

    side = tip_sets(true_tree)[true_tree.root.clades[0]]  # one side of the true root split
    m['root_recovered'] = any(sets[c] in (side, sets[tree.root] - side) for c in tree.root.clades)
    m['root_date_error'] = numdate[tree.root] - dates.loc['node_root', 'date']

    matched = [(c, true_by_set[sets[c]]) for c in tree.get_nonterminals() if sets[c] in true_by_set]
    m['frac_nodes_matched'] = len(matched) / len(tree.get_nonterminals())
    err = np.array([numdate[c] - dates.loc[n, 'date'] for c, n in matched])
    tau = np.array([dates.loc[n, 'tau'] for _, n in matched])
    m |= {'date_bias': float(err.mean()), 'date_rmse': float(np.sqrt(np.mean(err**2))), 'date_mae': float(np.abs(err).mean())}
    m |= binned(err, tau, 'date_rmse')

    ci = [(nodes[c.name].get('num_date_confidence'), dates.loc[n, 'date']) for c, n in matched]
    if all(x is not None for x, _ in ci):
        m['ci90_coverage'] = float(np.mean([lo <= t <= hi for (lo, hi), t in ci]))
        m['ci90_mean_width'] = float(np.mean([hi - lo for (lo, hi), _ in ci]))

    coal_tsv = out / 'coalescent.tsv'
    if coal_tsv.exists():
        seg = pd.read_csv(coal_tsv, sep='\t')
        seg_tau = PRESENT - 0.5 * (seg['segment.start'] + seg['segment.end'])
        seg = seg[(seg_tau >= 0) & (seg_tau <= dates.loc['node_root', 'tau'])]
        tc_true = SCHEDULES[coal]().tc(PRESENT - 0.5 * (seg['segment.start'] + seg['segment.end']))
        log_err = np.log(seg['T_c.value'] / tc_true)
        m['tc_log_bias'] = float(log_err.mean())
        m['tc_log_mae'] = float(np.abs(log_err).mean())
        m['tc_coverage'] = float(np.mean((seg['T_c.lower'] <= tc_true) & (tc_true <= seg['T_c.upper'])))
    return m


def evaluate(out, suite, data):
    cache = out / 'metrics.json'
    if cache.exists() and cache.stat().st_mtime >= (out / 'run.json').stat().st_mtime:
        return json.loads(cache.read_text())
    run = json.loads((out / 'run.json').read_text())
    variant, coal, rep, setname = out.parts[-4:]
    L, mu = setname.removeprefix('L').split('_mu')
    row = {'variant': variant, 'coal': coal, 'rep': rep, 'L': int(L), 'mu': float(mu)}
    row |= {'returncode': run['returncode'], 'seconds': run['seconds']}
    if run['returncode'] == 0:
        set_dir = data / coal / rep / setname
        row |= eval_ancestral(out, set_dir) if suite == 'ancestral' else eval_timetree(out, set_dir, coal)
    cache.write_text(json.dumps(row, indent=1, default=float))
    return row


SUMMARY = {
    'ancestral': ['site_error', 'pi_max_abs_error', 'W_log_rms_error', 'mu_inferred'],
    'timetree': ['clock_rate_rel_error', 'root_date_error', 'date_bias', 'date_rmse', 'n_negative_branches',
                 'root_recovered', 'frac_nodes_matched', 'tc_log_bias', 'tc_log_mae', 'ci90_coverage'],
}  # fmt: skip


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--suite', choices=['ancestral', 'timetree'], required=True)
    p.add_argument('--data', type=Path, default=HERE / 'data')
    p.add_argument('--results', type=Path, default=HERE / 'results')
    p.add_argument('-j', '--jobs', type=int, default=available_cpus())
    args = p.parse_args()

    outs = sorted(r.parent for r in (args.results / args.suite).glob('*/*/rep*/L*/run.json'))
    with ProcessPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(evaluate, outs, [args.suite] * len(outs), [args.data] * len(outs)))
    df = pd.DataFrame(rows)
    df.to_csv(args.results / f'{args.suite}.tsv', sep='\t', index=False)

    failed = df[df.returncode != 0]
    print(f'{len(df)} runs, {len(failed)} failed')
    cols = [c for c in SUMMARY[args.suite] if c in df]
    with pd.option_context('display.width', 250, 'display.max_columns', 30, 'display.float_format', '{:.4g}'.format):
        print(df[df.returncode == 0].groupby(['variant', 'coal', 'L', 'mu'])[cols].mean())


if __name__ == '__main__':
    main()
