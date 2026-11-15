"""Validation of TreeTime on segmented genomes with reassortment: simulate, run, evaluate.

  python segments.py simulate [--reps 3]      writes data_segments/<case>/rep<NN>/
  python segments.py run [-j 8]               runs ../.out/treetime, writes results_segments/<case>/rep<NN>/<variant>/
  python segments.py evaluate                 writes results_segments/branches.tsv, rates.tsv and summary.md

Data: an ancestral reassortment graph (reassortment.py, Tc = 10 years, 80 tips at 40 sample times, rate `rho` per
lineage and year) and four segments with their own length, substitution rate and random GTR. Sequences evolve along
each segment's local tree. TreeTime sees, per segment, the tip alignment and a divergence tree on the true topology
with branch lengths mu t * lognormal(0, 0.3), named internal nodes (perfect node identities, as from an ideal
TreeKnit), in two forms: `full` keeps the nodes that the segment only passes through, `collapsed` drops them, as tree
inference would.

Variants: `joint-full` and `joint-collapsed` run `optimize --segment-trees` on all segments, `separate` runs
`optimize --tree` once per segment (collapsed tree, topology kept).

Expected: the relative segment rates are recovered; on branches shared by several segments the joint fit is closer
to the true substitutions than the separate fits, the more so the more segments share the branch; branches of one
segment are fitted about as well as by the separate fit; high reassortment leaves fewer shared branches.
"""

import argparse
import json
import lzma
import subprocess
import sys
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

from coalescent import Constant
from reassortment import collapsed, local_tree, simulate_arg

HERE = Path(__file__).resolve().parent
SEGMENTS = {'s1': (1500, 4e-3), 's2': (1200, 3e-3), 's3': (1000, 2e-3), 's4': (800, 1e-3)}  # L, mu per site and year
CASES = {'low': 0.02, 'high': 0.2}  # reassortment rate per lineage and year
PERTURB_SIGMA = 0.3
KEEP_TOPOLOGY = ['--no-collapse-short-branches', '--no-merge-shared-mutations', '--no-hoist-reversions']
ALPHABET = np.array(list('ACGT'))


def random_gtr(rng):
    w = rng.random((4, 4))
    w = w + w.T
    np.fill_diagonal(w, 0)
    pi = rng.random(4)
    pi /= pi.sum()
    q = w * pi[None, :]
    np.fill_diagonal(q, -q.sum(axis=1))
    q /= -np.dot(pi, np.diag(q))
    return pi, w, q


def transition(q, pi, t):
    d = np.sqrt(pi)
    sym = d[:, None] * q / d[None, :]
    vals, vecs = np.linalg.eigh(sym)
    p = (vecs * np.exp(vals * t)) @ vecs.T
    p = np.clip(p * d[None, :] / d[:, None], 0, None)
    return p / p.sum(axis=1, keepdims=True)


def newick(root, children, length, name=lambda n: n):
    def rec(node):
        kids = children.get(node, [])
        inner = '(' + ','.join(rec(c) for c in kids) + ')' if kids else ''
        tail = '' if node == root else f':{length(node):.10g}'
        return f'{inner}{name(node)}{tail}'

    return rec(root) + ';\n'


def simulate(args):
    for case, rho in CASES.items():
        for rep in range(args.reps):
            out = args.data / case / f'rep{rep:02d}'
            if (out / 'truth.json').exists():
                continue
            rng = np.random.default_rng([args.seed, list(CASES).index(case), rep])
            arg = simulate_arg(Constant(10.0), rng, list(SEGMENTS), rho)
            out.mkdir(parents=True, exist_ok=True)
            truth = {'rho': rho, 'n_reassortments': arg.n_reassortments, 'segments': {}}
            for segment, (L, mu) in SEGMENTS.items():
                root, children, parent = local_tree(arg, segment)
                pi, w, q = random_gtr(rng)
                seqs = {root: rng.choice(4, size=L, p=pi)}
                stack = [root]
                while stack:
                    node = stack.pop()
                    for c in children.get(node, []):
                        p = transition(q, pi, mu * (arg.tau[node] - arg.tau[c]))
                        u = rng.random(L)[:, None]
                        seqs[c] = np.minimum((u > np.cumsum(p[seqs[node]], axis=1)).sum(axis=1), 3)
                        stack.append(c)
                with lzma.open(out / f'{segment}.tips.fasta.xz', 'wt') as fh:
                    for tip in arg.tips:
                        fh.write(f'>{tip}\n{"".join(ALPHABET[seqs[tip]])}\n')
                with lzma.open(out / f'{segment}.nodes.fasta.xz', 'wt') as fh:
                    for node, seq in seqs.items():
                        fh.write(f'>{node}\n{"".join(ALPHABET[seq])}\n')
                (out / f'{segment}.time.nwk').write_text(
                    newick(root, children, lambda n: arg.tau[parent[n]] - arg.tau[n])
                )
                perturb = {n: rng.lognormal(0.0, PERTURB_SIGMA) for n in parent}
                (out / f'{segment}.full.nwk').write_text(
                    newick(root, children, lambda n: mu * (arg.tau[parent[n]] - arg.tau[n]) * perturb[n])
                )
                croot, cchildren, cparent = collapsed(root, children)
                factor = {n: rng.lognormal(0.0, PERTURB_SIGMA) for n in cparent}
                (out / f'{segment}.collapsed.nwk').write_text(
                    newick(croot, cchildren, lambda n: mu * (arg.tau[cparent[n]] - arg.tau[n]) * factor[n])
                )
                truth['segments'][segment] = {'L': L, 'mu': mu, 'pi': pi.tolist(), 'W': w.tolist()}
            (out / 'truth.json').write_text(json.dumps(truth, indent=1))
            print(f'{out}: {arg.n_reassortments} reassortments', flush=True)


def jobs(args):
    for rep_dir in sorted(args.data.glob('*/rep*')):
        if not (rep_dir / 'truth.json').exists():
            continue
        out = args.results / rep_dir.relative_to(args.data)
        alignments = [str(rep_dir / f'{s}.tips.fasta.xz') for s in SEGMENTS]
        for form in ('full', 'collapsed'):
            trees = [str(rep_dir / f'{s}.{form}.nwk') for s in SEGMENTS]
            target = out / f'joint-{form}'
            yield target, [
                '--segment-trees', *trees, '--segment-alignments', *alignments,
                '--output-tree-nwk', str(target / 'tree.{alignment}.nwk'), '--output-model', str(target / 'model.json'),
            ]  # fmt: skip
        for index, segment in enumerate(SEGMENTS, 1):
            target = out / 'separate' / segment
            yield target, [
                '--tree', str(rep_dir / f'{segment}.collapsed.nwk'), '--alignment', alignments[index - 1], *KEEP_TOPOLOGY,
                '--output-tree-nwk', str(target / 'tree.nwk'), '--output-model', str(target / 'model.json'),
            ]  # fmt: skip


def run_one(binary, target, flags):
    if (target / 'done').exists():
        return 0
    target.mkdir(parents=True, exist_ok=True)
    with open(target / 'log.txt', 'w') as log:
        rc = subprocess.call([binary, 'optimize', '-v', *flags], stdout=log, stderr=subprocess.STDOUT)
    if rc == 0:
        (target / 'done').touch()
    print(f'{"ok " if rc == 0 else "FAIL"} {target}', flush=True)
    return rc


def run(args):
    with ThreadPoolExecutor(args.jobs) as pool:
        failed = sum(pool.map(lambda job: run_one(args.binary, *job), jobs(args)))
    print(f'{failed} failed')


def read_fasta_xz(path):
    seqs, name = {}, None
    with lzma.open(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:]
                seqs[name] = []
            else:
                seqs[name].append(line)
    return {k: np.frombuffer(''.join(v).encode(), 'S1') for k, v in seqs.items()}


def parents_with_lengths(path):
    tree = Phylo.read(path, 'newick')
    return {c.name: (p.name, c.branch_length or 0.0) for p in tree.find_clades() for c in p.clades}


def path_length(edges, child, ancestor):
    total, node = 0.0, child
    while node != ancestor:
        node, length = edges[node]
        total += length
    return total


def model_rates(path):
    return [m['mu'] for m in json.loads(path.read_text())['models']]


def evaluate(args):
    rows, rate_rows = [], []
    for rep_dir in sorted(args.data.glob('*/rep*')):
        truth = json.loads((rep_dir / 'truth.json').read_text())
        out = args.results / rep_dir.relative_to(args.data)
        collapsed_edges = {s: parents_with_lengths(rep_dir / f'{s}.collapsed.nwk') for s in SEGMENTS}
        time_edges = {s: parents_with_lengths(rep_dir / f'{s}.time.nwk') for s in SEGMENTS}
        shared = {}
        for s, edges in collapsed_edges.items():
            for child, (parent, _) in edges.items():
                shared[(parent, child)] = shared.get((parent, child), 0) + 1
        fitted = {}
        for variant in ('joint-full', 'joint-collapsed'):
            rates = model_rates(out / variant / 'model.json')
            fitted[variant] = {
                s: (parents_with_lengths(out / variant / f'tree.{i}.nwk'), rates[i - 1])
                for i, s in enumerate(SEGMENTS, 1)
            }
        fitted['separate'] = {
            s: (parents_with_lengths(out / 'separate' / s / 'tree.nwk'), model_rates(out / 'separate' / s / 'model.json')[0])
            for s in SEGMENTS
        }
        lengths = np.array([SEGMENTS[s][0] for s in SEGMENTS], dtype=float)
        true_rates = np.array([SEGMENTS[s][1] for s in SEGMENTS])
        for variant in ('joint-full', 'joint-collapsed'):
            est = np.array([fitted[variant][s][1] for s in SEGMENTS])
            for s, t, e in zip(SEGMENTS, true_rates / np.average(true_rates, weights=lengths), est / np.average(est, weights=lengths)):
                rate_rows.append({'case': rep_dir.parent.name, 'rep': rep_dir.name, 'variant': variant, 'segment': s,
                                  'true_relative_rate': t, 'estimated_relative_rate': e})  # fmt: skip
        for s, (L, mu) in SEGMENTS.items():
            seqs = read_fasta_xz(rep_dir / f'{s}.nodes.fasta.xz')
            parents = {parent for parent, _ in collapsed_edges[s].values()}
            root = next(p for p in parents if p not in collapsed_edges[s])
            for child, (parent, _) in collapsed_edges[s].items():
                time = path_length(time_edges[s], child, parent)
                realized = np.count_nonzero(seqs[child] != seqs[parent])
                row = {'case': rep_dir.parent.name, 'rep': rep_dir.name, 'segment': s, 'parent': parent, 'child': child,
                       'n_segments': shared[(parent, child)], 'at_root': parent == root,
                       'expected': mu * time * L, 'realized': realized}  # fmt: skip
                for variant, per_segment in fitted.items():
                    edges, rate = per_segment[s]
                    row[variant] = rate * path_length(edges, child, parent) * L
                rows.append(row)
    branches = pd.DataFrame(rows)
    rates = pd.DataFrame(rate_rows)
    args.results.mkdir(parents=True, exist_ok=True)
    branches.to_csv(args.results / 'branches.tsv', sep='\t', index=False, float_format='%.6g')
    rates.to_csv(args.results / 'rates.tsv', sep='\t', index=False, float_format='%.6g')
    summary(branches, rates, args.data, args.results / 'summary.md')


def summary(branches, rates, data, path):
    lines = ['# Segment validation', '']
    meta = [json.loads(p.read_text()) | {'case': p.parent.parent.name} for p in sorted(data.glob('*/rep*/truth.json'))]
    lines += ['Reassortments per replicate: ' + ', '.join(f'{m["case"]} {m["n_reassortments"]}' for m in meta), '']
    lines += ['## Relative segment rates (length-weighted mean 1)', '']
    table = rates.groupby(['variant', 'segment']).agg(
        true=('true_relative_rate', 'mean'), estimated=('estimated_relative_rate', 'mean'),
        sd=('estimated_relative_rate', 'std'),
    )  # fmt: skip
    lines += [table.round(3).to_markdown(), '']
    variants = ['separate', 'joint-collapsed', 'joint-full']
    for reference, note in [
        ('expected', 'the expected substitutions mu t L, which a fit pooled over segments targets'),
        ('realized', "the segment's own realized substitutions, which a separate fit follows"),
    ]:
        lines += [f'## Branch error against {reference} substitutions', '',
                  f'RMSE in substitutions of the fitted substitutions (rate x length x L) against {note}, over the '
                  'branches of the collapsed segment trees, by the number of segments whose collapsed tree contains '
                  'the same branch. `poisson` is the RMSE of Poisson counts around mu t L.', '']  # fmt: skip
        grouped = branches[~branches['at_root']].groupby(['case', 'n_segments'])
        table = pd.DataFrame({'branches': grouped.size()})
        for variant in variants:
            table[variant] = grouped.apply(lambda g, v=variant: np.sqrt(np.mean((g[v] - g[reference]) ** 2)))
        table['poisson'] = np.sqrt(grouped['expected'].mean())
        lines += [table.round(2).to_markdown(), '']
    roots = branches[branches['at_root']].groupby(['case', 'rep', 'segment'])[['expected', 'realized', *variants]].sum()
    roots = roots.groupby('case').apply(
        lambda g: pd.Series({v: np.sqrt(np.mean((g[v] - g['expected']) ** 2)) for v in variants})
    )
    lines += ['## Sum of the two branches at each segment root against expected substitutions (RMSE)', '',
              roots.round(2).to_markdown(), '']  # fmt: skip
    path.write_text('\n'.join(lines))
    print('\n'.join(lines))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('step', choices=['simulate', 'run', 'evaluate'])
    p.add_argument('--data', type=Path, default=HERE / 'data_segments')
    p.add_argument('--results', type=Path, default=HERE / 'results_segments')
    p.add_argument('--reps', type=int, default=3)
    p.add_argument('--seed', type=int, default=2026)
    p.add_argument('--binary', default=str(HERE.parent / '.out' / 'treetime'))
    p.add_argument('-j', '--jobs', type=int, default=4)
    args = p.parse_args()
    {'simulate': simulate, 'run': run, 'evaluate': evaluate}[args.step](args)


if __name__ == '__main__':
    sys.exit(main())
