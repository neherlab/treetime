"""Validation of TreeTime on segmented genomes with reassortment: simulate, run, evaluate.

  python segments.py simulate [--reps 3]      writes data_segments/<case>/rep<NN>/
  python segments.py run [-j 8]               runs ../.out/treetime, writes results_segments/<case>/rep<NN>/<variant>/
  python segments.py evaluate                 writes results_segments/branches.tsv, rates.tsv and summary.md
  python segments.py export --case high --rep 0 --out ../data/simulated/reassortment/80
                                              writes one replicate as an example dataset of TreeTime

Data: an ancestral reassortment graph (reassortment.py, Tc = 10 years, 80 tips at 40 sample times, rate `rho` per
lineage and year) and four segments with their own length, substitution rate and random GTR. Sequences evolve along
each segment's local tree. TreeTime sees, per segment, the tip alignment and a divergence tree on the true topology
with branch lengths mu t * lognormal(0, 0.3), named internal nodes (perfect node identities, as from an ideal
TreeKnit), in three forms: `full` keeps the nodes that the segment only passes through where another segment branches,
`collapsed` drops them, as tree inference would. These two forms leave out the hidden nodes, which have one child in
every segment (reassortment nodes, and mergers of lineages with disjoint segments): no tree inference recovers them,
and an edge across a dropped hidden node gets the sum of the perturbed lengths. `kept` is `full` with the hidden nodes,
as nodes with one child in every segment, and the same perturbed lengths; it shows what hidden nodes in the input do.

Variants: `joint-kept`, `joint-full` and `joint-collapsed` run `optimize --segment-trees` on all segments, `separate`
runs `optimize --tree` once per segment (collapsed tree, topology kept). `timetree-kept`, `timetree-full` and
`timetree-collapsed` run `timetree --segment-trees` at the true mean clock rate with the default message passing, where
the message to a parent is the subtree alone, the `-cross` variants with `--cross-parent-messages`. Each timetree run
also writes the merged graph with its fitted divergence lengths (`graph.dot`), which the evaluation compares with the
input trees.

Expected: the relative segment rates are recovered; on branches shared by several segments the joint fit is closer
to the true substitutions than the separate fits, the more so the more segments share the branch; branches of one
segment are fitted about as well as by the separate fit; high reassortment leaves fewer shared branches. The dates of internal nodes are unbiased, and their 90% regions
cover about 90% of the true dates; cross-parent messages matter at nodes with several parents.
"""

import argparse
import json
import lzma
import re
import shutil
import subprocess
import sys
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

from coalescent import Constant
from coalescent import PRESENT
from reassortment import collapsed, hidden_nodes, local_tree, simulate_arg

HERE = Path(__file__).resolve().parent
SEGMENTS = {'s1': (1500, 4e-3), 's2': (1200, 3e-3), 's3': (1000, 2e-3), 's4': (800, 1e-3)}  # L, mu per site and year
CASES = {'low': 0.02, 'high': 0.2}  # reassortment rate per lineage and year
PERTURB_SIGMA = 0.3
KEEP_TOPOLOGY = ['--no-collapse-short-branches', '--no-merge-shared-mutations', '--no-hoist-reversions']
ALPHABET = np.array(list('ACGT'))
SAMPLES_PER_BATCH, SAMPLE_INTERVAL = 2, 0.1  # the sampling schedule of simulate_arg
TREE_FORMS = ('kept', 'full', 'collapsed')
TIMETREE_VARIANTS = {
    'timetree-kept': ('kept', []),
    'timetree-full': ('full', []),
    'timetree-full-cross': ('full', ['--cross-parent-messages']),
    'timetree-collapsed': ('collapsed', []),
    'timetree-collapsed-cross': ('collapsed', ['--cross-parent-messages']),
}


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
            if all((out / name).exists() for name in ('truth.json', *(f'{s}.kept.nwk' for s in SEGMENTS))):
                continue
            rng = np.random.default_rng([args.seed, list(CASES).index(case), rep])
            arg = simulate_arg(Constant(10.0), rng, list(SEGMENTS), rho)
            out.mkdir(parents=True, exist_ok=True)
            truth = {'rho': rho, 'n_reassortments': arg.n_reassortments, 'segments': {}}
            local_trees = {segment: local_tree(arg, segment) for segment in SEGMENTS}
            hidden = hidden_nodes(local_trees.values())
            for segment, (L, mu) in SEGMENTS.items():
                root, children, parent = local_trees[segment]
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
                (out / f'{segment}.kept.nwk').write_text(
                    newick(root, children, lambda n: mu * (arg.tau[parent[n]] - arg.tau[n]) * perturb[n])
                )
                froot, fchildren, fparent = collapsed(root, children, hidden)
                edges = {n: (p, mu * (arg.tau[p] - arg.tau[n]) * perturb[n]) for n, p in parent.items()}
                (out / f'{segment}.full.nwk').write_text(
                    newick(froot, fchildren, lambda n: path_length(edges, n, fparent[n]))
                )
                croot, cchildren, cparent = collapsed(root, children)
                factor = {n: rng.lognormal(0.0, PERTURB_SIGMA) for n in cparent}
                (out / f'{segment}.collapsed.nwk').write_text(
                    newick(croot, cchildren, lambda n: mu * (arg.tau[cparent[n]] - arg.tau[n]) * factor[n])
                )
                truth['segments'][segment] = {'L': L, 'mu': mu, 'pi': pi.tolist(), 'W': w.tolist()}
            (out / 'truth.json').write_text(json.dumps(truth, indent=1))
            print(f'{out}: {arg.n_reassortments} reassortments', flush=True)


def mean_clock_rate():
    lengths = np.array([L for L, _ in SEGMENTS.values()], dtype=float)
    return float(np.average([mu for _, mu in SEGMENTS.values()], weights=lengths))


def tip_date(tip):
    index = int(tip.split('_')[1]) - 1
    return PRESENT - (index // SAMPLES_PER_BATCH) * SAMPLE_INTERVAL


def node_dates(rep_dir):
    dates = {}
    for s in SEGMENTS:
        tree = Phylo.read(rep_dir / f'{s}.time.nwk', 'newick')
        depths = tree.depths()
        tip = next(t for t in tree.get_terminals())
        anchor = tip_date(tip.name) - depths[tip]
        dates.update({clade.name: anchor + depth for clade, depth in depths.items()})
    return dates


def write_metadata(rep_dir):
    path = rep_dir / 'metadata.tsv'
    if not path.exists():
        tips = sorted(name for name in node_dates(rep_dir) if name.startswith('tip_'))
        path.write_text('name\tdate\n' + ''.join(f'{tip}\t{tip_date(tip):.4f}\n' for tip in tips))
    return path


def jobs(args):
    for rep_dir in sorted(args.data.glob('*/rep*')):
        if not (rep_dir / 'truth.json').exists():
            continue
        out = args.results / rep_dir.relative_to(args.data)
        alignments = [str(rep_dir / f'{s}.tips.fasta.xz') for s in SEGMENTS]
        metadata = write_metadata(rep_dir)
        for variant, (form, extra) in TIMETREE_VARIANTS.items():
            trees = [str(rep_dir / f'{s}.{form}.nwk') for s in SEGMENTS]
            target = out / variant
            yield target, 'timetree', [
                '--segment-trees', *trees, '--segment-alignments', *alignments, '--metadata', str(metadata),
                '--clock-rate', f'{mean_clock_rate():.6g}', '--confidence', *extra,
                '--output-dates', str(target / 'dates.tsv'), '--output-tree-dot', str(target / 'graph.dot'),
            ]  # fmt: skip
        for form in TREE_FORMS:
            trees = [str(rep_dir / f'{s}.{form}.nwk') for s in SEGMENTS]
            target = out / f'joint-{form}'
            yield target, 'optimize', [
                '--segment-trees', *trees, '--segment-alignments', *alignments,
                '--output-tree-nwk', str(target / 'tree.{alignment}.nwk'), '--output-model', str(target / 'model.json'),
            ]  # fmt: skip
        for index, segment in enumerate(SEGMENTS, 1):
            target = out / 'separate' / segment
            yield target, 'optimize', [
                '--tree', str(rep_dir / f'{segment}.collapsed.nwk'), '--alignment', alignments[index - 1], *KEEP_TOPOLOGY,
                '--output-tree-nwk', str(target / 'tree.nwk'), '--output-model', str(target / 'model.json'),
            ]  # fmt: skip


def run_one(binary, target, command, flags):
    if (target / 'done').exists():
        return 0
    target.mkdir(parents=True, exist_ok=True)
    with open(target / 'log.txt', 'w') as log:
        rc = subprocess.call([binary, command, '-v', *flags], stdout=log, stderr=subprocess.STDOUT)
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
        for variant in (f'joint-{form}' for form in TREE_FORMS):
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
        for variant in (f'joint-{form}' for form in TREE_FORMS):
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
    dates = evaluate_dates(args)
    dates.to_csv(args.results / 'dates.tsv', sep='\t', index=False, float_format='%.6g')
    hidden_edges = evaluate_hidden_edges(args)
    hidden_edges.to_csv(args.results / 'hidden_edges.tsv', sep='\t', index=False, float_format='%.6g')
    hidden_paths = evaluate_hidden_paths(args)
    hidden_paths.to_csv(args.results / 'hidden_paths.tsv', sep='\t', index=False, float_format='%.6g')
    graphs = evaluate_graphs(args)
    graphs.to_csv(args.results / 'graphs.tsv', sep='\t', index=False)
    hidden = {'edges': hidden_edges, 'paths': hidden_paths, 'graphs': graphs}
    summary(branches, rates, dates, hidden, args.data, args.results / 'summary.md')


def read_dates(path):
    table = pd.read_csv(path, sep='\t', comment=None, skiprows=1)
    table.columns = ['node', 'date', 'numeric', 'lower', 'upper']
    return table.set_index('node')


def node_kinds(rep_dir, form):
    """'hidden' for nodes with one child in every segment, else 'several parents' or 'one parent' in the graph."""
    parents, branching = {}, set()
    for s in SEGMENTS:
        edges = parents_with_lengths(rep_dir / f'{s}.{form}.nwk')
        for child, (parent, _) in edges.items():
            parents.setdefault(child, set()).add(parent)
        counts = pd.Series([parent for parent, _ in edges.values()]).value_counts()
        branching |= set(counts[counts > 1].index)
    return {
        node: 'hidden' if node not in branching else 'several parents' if len(p) > 1 else 'one parent'
        for node, p in parents.items()
    } | {node: 'one parent' for node in branching if node not in parents}


def graph_edges(rep_dir, form):
    """{(parent, child): segments} of the graph that the segment trees of `form` give when merged by name."""
    edges = {}
    for s in SEGMENTS:
        for child, (parent, _) in parents_with_lengths(rep_dir / f'{s}.{form}.nwk').items():
            edges.setdefault((parent, child), set()).add(s)
    return edges


def hidden_kinds(rep_dir):
    """Kind of each hidden node of the `kept` trees: 'reassortment' (several parents), 'disjoint merger' (several
    children), or 'loop' (one parent and one child: a reassortment whose two lineages merge again, and that merger)."""
    edges = graph_edges(rep_dir, 'kept')
    parents, children = {}, {}
    for parent, child in edges:
        parents.setdefault(child, set()).add(parent)
        children.setdefault(parent, set()).add(child)
    hidden = [node for node, kind in node_kinds(rep_dir, 'kept').items() if kind == 'hidden' and node in children]
    return {
        node: 'reassortment' if len(parents[node]) > 1 else 'disjoint merger' if len(children[node]) > 1 else 'loop'
        for node in hidden
    }


def hidden_neighbours(rep_dir, hidden):
    """Nodes that are not hidden and are linked to a hidden node by an edge or by a path through hidden nodes."""
    adjacent = {}
    for parent, child in graph_edges(rep_dir, 'kept'):
        adjacent.setdefault(parent, set()).add(child)
        adjacent.setdefault(child, set()).add(parent)
    reached = set().union(*(adjacent[node] for node in hidden)) if hidden else set()
    return reached - set(hidden)


def evaluate_dates(args):
    rows = []
    for rep_dir in sorted(args.data.glob('*/rep*')):
        truth = node_dates(rep_dir)
        kinds = {form: node_kinds(rep_dir, form) for form in TREE_FORMS}
        hidden = hidden_kinds(rep_dir)
        neighbours = hidden_neighbours(rep_dir, hidden)
        out = args.results / rep_dir.relative_to(args.data)
        for variant in TIMETREE_VARIANTS:
            path = out / variant / 'dates.tsv'
            if not path.exists():
                continue
            log = (out / variant / 'log.txt').read_text()
            form = TIMETREE_VARIANTS[variant][0]
            sweeps = next((int(line.split(' took ')[1].split()[0]) for line in log.splitlines() if ' took ' in line and 'sweeps' in line), 1)
            inferred = read_dates(path)
            for node, date in truth.items():
                if node.startswith('tip_') or node not in inferred.index:
                    continue
                row = inferred.loc[node]
                rows.append({'case': rep_dir.parent.name, 'rep': rep_dir.name, 'variant': variant, 'node': node,
                             'kind': kinds[form].get(node, 'hidden'), 'full_kind': kinds['full'].get(node, ''),
                             'hidden_kind': hidden.get(node, ''), 'neighbour': node in neighbours,
                             'true': date, 'inferred': float(row['numeric']),
                             'lower': float(row['lower']), 'upper': float(row['upper']), 'sweeps': sweeps})  # fmt: skip
    return pd.DataFrame(rows)


def read_dot(path):
    """{(parent, child): divergence} of a DOT graph that TreeTime writes."""
    names, lengths = {}, {}
    for line in path.read_text().splitlines():
        if (m := re.match(r'\s*(\d+) \[label="\(\d+\) (\S+)"\]', line)) is not None:
            names[m[1]] = m[2]
        elif 'invis' not in line and (m := re.match(r'\s*(\d+) -> (\d+) \[xlabel=.*weight="([^"]+)"\]', line)):
            lengths[(m[1], m[2])] = float(m[3])
    return {(names[p], names[c]): length for (p, c), length in lengths.items()}


def contract(edges, hidden):
    """`edges` ({(parent, child): segments}) without the `hidden` nodes: each child takes over the parent edges whose
    segments it shares."""
    edges = dict(edges)
    for node in hidden:
        parents = [(p, segs) for (p, c), segs in edges.items() if c == node]
        children = [(c, segs) for (p, c), segs in edges.items() if p == node]
        edges = {key: segs for key, segs in edges.items() if node not in key}
        for p, parent_segments in parents:
            for c, child_segments in children:
                if common := parent_segments & child_segments:
                    edges[(p, c)] = edges.get((p, c), set()) | common
    return edges


def evaluate_graphs(args):
    rows = []
    for rep_dir in sorted(args.data.glob('*/rep*')):
        out = args.results / rep_dir.relative_to(args.data)
        hidden = hidden_kinds(rep_dir)
        given = {form: graph_edges(rep_dir, form) for form in ('kept', 'full')}
        dots = {form: out / f'timetree-{form}' / 'graph.dot' for form in ('kept', 'full')}
        merged = {form: read_dot(path) for form, path in dots.items() if path.exists()}
        row = {'case': rep_dir.parent.name, 'rep': rep_dir.name, 'hidden': len(hidden),
               'input_contraction_equal': contract(given['kept'], hidden) == given['full']}  # fmt: skip
        for form, graph in merged.items():
            row[f'{form}_missing'] = ' '.join(f'{p}->{c}' for p, c in sorted(set(given[form]) - set(graph))) or '-'
        if len(merged) == 2:
            contracted = set(contract({key: given['kept'][key] for key in merged['kept']}, hidden))
            row['merged_contraction_differs'] = len(contracted ^ set(merged['full']))
        if 'kept' in merged:
            below = [length for (p, _), length in merged['kept'].items() if p in hidden]
            row['below_hidden_zero'] = f'{sum(length < 1e-8 for length in below)}/{len(below)}'
        rows.append(row)
    return pd.DataFrame(rows)


def evaluate_hidden_paths(args):
    """For each parent and child joined through hidden nodes in the kept graph: the true, fitted and inferred spans
    in kept, and the fitted and inferred spans of the same pair in contracted."""
    rows = []
    rate = mean_clock_rate()
    for rep_dir in sorted(args.data.glob('*/rep*')):
        out = args.results / rep_dir.relative_to(args.data)
        variants = {form: out / f'timetree-{form}' for form in ('kept', 'full')}
        if not all((v / 'graph.dot').exists() for v in variants.values()):
            continue
        truth, hidden = node_dates(rep_dir), hidden_kinds(rep_dir)
        fitted = {form: read_dot(v / 'graph.dot') for form, v in variants.items()}
        dates = {form: read_dates(v / 'dates.tsv')['numeric'] for form, v in variants.items()}
        for path in hidden_paths(graph_edges(rep_dir, 'kept'), hidden):
            parent, child = path[0], path[-1]
            kept_divergence = sum(fitted['kept'].get(edge, np.nan) for edge in zip(path, path[1:]))
            rows.append({
                'case': rep_dir.parent.name, 'rep': rep_dir.name, 'parent': parent, 'child': child,
                'hidden_on_path': len(path) - 2, 'true_span': truth[child] - truth[parent],
                'kept_fitted_span': kept_divergence / rate, 'kept_span': dates['kept'][child] - dates['kept'][parent],
                'contracted_fitted_span': fitted['full'].get((parent, child), np.nan) / rate,
                'contracted_span': dates['full'][child] - dates['full'][parent],
            })  # fmt: skip
    return pd.DataFrame(rows).drop_duplicates(['case', 'rep', 'parent', 'child']) if rows else pd.DataFrame()


def hidden_paths(edges, hidden):
    """Paths [parent, hidden..., child] whose inner nodes are hidden and whose edges share a segment."""
    children = {}
    for (p, c), segments in edges.items():
        children.setdefault(p, []).append((c, segments))

    def walk(path, segments):
        for c, s in children.get(path[-1], []):
            if common := s & segments:
                yield from walk([*path, c], common) if c in hidden else [[*path, c]]

    for parent in sorted({p for (p, c) in edges if c in hidden and p not in hidden}):
        for c, s in children[parent]:
            if c in hidden:
                yield from walk([parent, c], s)


def evaluate_hidden_edges(args):
    """Substitutions that `joint-kept` fits on the edges below and above each hidden node, per segment of the edge,
    against the expected substitutions mu t L."""
    rows = []
    for rep_dir in sorted(args.data.glob('*/rep*')):
        fitted_dir = args.results / rep_dir.relative_to(args.data) / 'joint-kept'
        if not (fitted_dir / 'model.json').exists():
            continue
        truth = node_dates(rep_dir)
        hidden = hidden_kinds(rep_dir)
        rates = model_rates(fitted_dir / 'model.json')
        fitted = {s: parents_with_lengths(fitted_dir / f'tree.{i}.nwk') for i, s in enumerate(SEGMENTS, 1)}
        for (parent, child), segments in graph_edges(rep_dir, 'kept').items():
            role = 'below' if parent in hidden else 'above' if child in hidden else None
            if role is None:
                continue
            node = parent if role == 'below' else child
            for i, s in enumerate(SEGMENTS):
                if s not in segments:
                    continue
                L, mu = SEGMENTS[s]
                fitted_parent, length = fitted[s].get(child, (None, np.nan))
                rows.append({'case': rep_dir.parent.name, 'rep': rep_dir.name, 'hidden': node, 'hidden_kind': hidden[node],
                             'role': role, 'parent': parent, 'child': child, 'segment': s,
                             'expected': mu * (truth[child] - truth[parent]) * L,
                             'fitted': rates[i] * length * L if fitted_parent == parent else np.nan})  # fmt: skip
    return pd.DataFrame(rows)


def summary(branches, rates, dates, hidden, data, path):
    lines = ['# Segment validation', '']
    meta = [json.loads(p.read_text()) | {'case': p.parent.parent.name} for p in sorted(data.glob('*/rep*/truth.json'))]
    lines += ['Reassortments per replicate: ' + ', '.join(f'{m["case"]} {m["n_reassortments"]}' for m in meta), '']
    lines += ['## Relative segment rates (length-weighted mean 1)', '']
    table = rates.groupby(['variant', 'segment']).agg(
        true=('true_relative_rate', 'mean'), estimated=('estimated_relative_rate', 'mean'),
        sd=('estimated_relative_rate', 'std'),
    )  # fmt: skip
    lines += [table.round(3).to_markdown(), '']
    variants = ['separate', 'joint-collapsed', 'joint-full', 'joint-kept']
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
    if not dates.empty:
        lines += dates_summary(dates)
        lines += hidden_summary(dates, hidden)
    path.write_text('\n'.join(lines))
    print('\n'.join(lines))


def dates_summary(dates):
    dates = dates.assign(
        error=dates['inferred'] - dates['true'],
        covered=(dates['lower'] <= dates['true']) & (dates['true'] <= dates['upper']),
        width=dates['upper'] - dates['lower'],
    )
    lines = ['## Dates of internal nodes', '',
             f'Fixed clock rate {mean_clock_rate():.4g} per site and year (the true length-weighted mean). Error is '
             'inferred minus true date in years; coverage is the fraction of true dates inside the 90% region. '
             '`sweeps` is the mean number of message-passing sweeps of the first round. Nodes with several parents are '
             'counted in the graph of the variant: the collapsed trees, which also drop the pass-through nodes, '
             'give more of them.', '']  # fmt: skip
    grouped = dates.groupby(['case', 'variant', 'kind'])
    table = pd.DataFrame({
        'nodes': grouped.size(),
        'bias': grouped['error'].mean(),
        'rmse': grouped['error'].apply(lambda e: np.sqrt(np.mean(e ** 2))),
        'coverage': grouped['covered'].mean(),
        'width': grouped['width'].median(),
        'sweeps': grouped['sweeps'].mean(),
    })  # fmt: skip
    return lines + [table.round(3).to_markdown(), '']


def hidden_summary(dates, hidden):
    forms = {'timetree-kept': 'kept', 'timetree-full': 'contracted', 'timetree-collapsed': 'collapsed'}
    dates = dates[dates['variant'].isin(forms)].assign(
        form=lambda d: d['variant'].map(forms),
        error=lambda d: d['inferred'] - d['true'],
        covered=lambda d: (d['lower'] <= d['true']) & (d['true'] <= d['upper']),
        width=lambda d: d['upper'] - d['lower'],
        age=lambda d: PRESENT - d['true'],
    )
    if dates['form'].nunique() < len(forms):
        return []
    counts = dates.groupby(['case', 'rep', 'node'])['form'].transform('nunique')
    shared = dates[(counts == len(forms)) & (dates['full_kind'] != '')]

    def stats(grouped):
        return pd.DataFrame({
            'nodes': grouped.size(),
            'bias': grouped['error'].mean(),
            'rmse': grouped['error'].apply(lambda e: np.sqrt(np.mean(e ** 2))),
            'coverage': grouped['covered'].mean(),
            'width': grouped['width'].median(),
        })  # fmt: skip

    lines = ['## Hidden nodes', '',
             '`kept` dates the `full` trees with the hidden nodes, `contracted` the `full` trees without them (the '
             '`timetree-full` variant), `collapsed` the collapsed trees. Shared nodes are the nodes that all three '
             'graphs date; the kind of a shared node is its kind in the contracted graph, so each row has the same '
             'nodes in every form.', '']  # fmt: skip
    table = stats(shared.groupby(['case', 'full_kind', 'form']))
    lines += ['### Shared nodes', '', table.round(3).to_markdown(), '']
    hidden_dates = dates[(dates['form'] == 'kept') & (dates['hidden_kind'] != '')]
    table = stats(hidden_dates.groupby(['case', 'hidden_kind']))
    lines += ['### Hidden nodes in `kept`', '', table.round(3).to_markdown(), '']
    pivot = shared.pivot_table(index=['case', 'rep', 'node'], columns='form', values='inferred')
    info = shared[shared['form'] == 'kept'].set_index(['case', 'rep', 'node'])[['true', 'age', 'neighbour', 'full_kind']]
    shift = info.join(pivot).assign(
        shift=lambda d: d['kept'] - d['contracted'],
        kept_error=lambda d: d['kept'] - d['true'],
        contracted_error=lambda d: d['contracted'] - d['true'],
        age_bin=lambda d: pd.cut(d['age'], [0, 5, 10, 20, np.inf], right=False),
    )
    grouped = shift.groupby(['case', 'neighbour', 'age_bin'], observed=True)
    table = pd.DataFrame({
        'nodes': grouped.size(),
        'rms_shift': grouped['shift'].apply(lambda e: np.sqrt(np.mean(e ** 2))),
        'max_abs_shift': grouped['shift'].apply(lambda e: e.abs().max()),
        'rmse_kept': grouped['kept_error'].apply(lambda e: np.sqrt(np.mean(e ** 2))),
        'rmse_contracted': grouped['contracted_error'].apply(lambda e: np.sqrt(np.mean(e ** 2))),
    })  # fmt: skip
    lines += ['### Date shift kept minus contracted at shared nodes', '',
              '`neighbour`: linked to a hidden node by an edge or through other hidden nodes. Age in years before '
              f'{PRESENT}.', '', table.round(3).to_markdown(), '']  # fmt: skip
    largest = shift.reindex(shift['shift'].abs().sort_values(ascending=False).index).head(12)
    columns = ['age', 'neighbour', 'full_kind', 'shift', 'kept_error', 'contracted_error']
    lines += ['Largest shifts:', '', largest[columns].round(3).to_markdown(), '']
    hidden_edges = hidden['edges']
    if not hidden_edges.empty:
        edges = hidden_edges.dropna(subset=['fitted'])
        grouped = edges.groupby(['case', 'role'])
        table = pd.DataFrame({
            'edges': grouped.size(),
            'fitted_below_0.01': grouped['fitted'].apply(lambda f: np.mean(f < 0.01)),
            'median_fitted': grouped['fitted'].median(),
            'median_expected': grouped['expected'].median(),
        })  # fmt: skip
        lines += ['### Edges at hidden nodes after `optimize --segment-trees` (`joint-kept`)', '',
                  'Fitted and expected substitutions per segment of the edge (rate x length x L, and mu t L). '
                  f'{len(hidden_edges) - len(edges)} edges are missing from the fitted trees.', '',
                  table.round(3).to_markdown(), '']  # fmt: skip
    paths = hidden['paths']
    if not paths.empty:
        paths = paths.assign(
            kept_stretch=paths['kept_span'] - paths['kept_fitted_span'],
            contracted_stretch=paths['contracted_span'] - paths['contracted_fitted_span'],
            kept_error=paths['kept_span'] - paths['true_span'],
            contracted_error=paths['contracted_span'] - paths['true_span'],
        )
        columns = ['kept_stretch', 'contracted_stretch', 'kept_error', 'contracted_error']
        table = paths.groupby(['case', 'hidden_on_path'])[columns].agg(['mean', 'median'])
        table.insert(0, 'paths', paths.groupby(['case', 'hidden_on_path']).size())
        lines += ['### Time spans of paths through hidden nodes', '',
                  'For each parent and child joined by a path through hidden nodes in the kept graph: stretch is the '
                  'inferred time span minus the fitted divergence of the path at the clock rate, in kept (sum over the '
                  'path) and in contracted (the edge from the parent to the child); error is the inferred minus the true '
                  'span. Years.', '', table.round(2).to_markdown(), '']  # fmt: skip
    graphs = hidden['graphs']
    if not graphs.empty:
        lines += ['### Merged graphs', '',
                  '`missing`: input edges that the merge of TreeTime rerouted; the `full` and `kept` trees hold every '
                  'pass-through node, so it should reroute none. `input_contraction_equal`: contracting the hidden nodes '
                  'of the kept trees, where each child takes over the parent edges whose segments it shares, gives the '
                  'edges and scopes of the full trees. `merged_contraction_differs`: edges in which the merged kept graph, '
                  'contracted, differs from the merged full graph. `below_hidden_zero`: edges below hidden nodes that the '
                  'optimize step of timetree-kept fits below 1e-8.', '',
                  graphs.to_markdown(index=False), '']  # fmt: skip
    return lines


def export(args):
    """One replicate as an example dataset: per segment the collapsed tree (`tree_<segment>.nwk`) and the tip alignment
    (`aln_<segment>.fasta.xz`), the tip dates (`metadata.tsv`), the true date of every node of the merged trees and its
    number of parents in the true graph without hidden nodes (`truth_nodes.tsv`), and the true model
    (`truth_model.json`)."""
    rep_dir = args.data / args.case / f'rep{args.rep:02d}'
    args.out.mkdir(parents=True, exist_ok=True)
    for s in SEGMENTS:
        shutil.copyfile(rep_dir / f'{s}.collapsed.nwk', args.out / f'tree_{s}.nwk')
        shutil.copyfile(rep_dir / f'{s}.tips.fasta.xz', args.out / f'aln_{s}.fasta.xz')
    shutil.copyfile(write_metadata(rep_dir), args.out / 'metadata.tsv')
    truth = json.loads((rep_dir / 'truth.json').read_text())
    model = {'clock_rate': mean_clock_rate(), **truth}
    (args.out / 'truth_model.json').write_text(json.dumps(model, indent=1) + '\n')
    parents = {}
    for parent, child in graph_edges(rep_dir, 'full'):
        parents.setdefault(child, set()).add(parent)
    nodes = sorted({node for edge in graph_edges(rep_dir, 'collapsed') for node in edge})
    dates = node_dates(rep_dir)
    rows = ''.join(f'{node}\t{dates[node]:.6f}\t{len(parents.get(node, ()))}\n' for node in nodes)
    (args.out / 'truth_nodes.tsv').write_text('name\tnumeric_date\tparents\n' + rows)
    print(f'{args.out}: {truth["n_reassortments"]} reassortments, {len(nodes)} nodes', flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('step', choices=['simulate', 'run', 'evaluate', 'export'])
    p.add_argument('--data', type=Path, default=HERE / 'data_segments')
    p.add_argument('--results', type=Path, default=HERE / 'results_segments')
    p.add_argument('--reps', type=int, default=3)
    p.add_argument('--seed', type=int, default=2026)
    p.add_argument('--binary', default=str(HERE.parent / '.out' / 'treetime'))
    p.add_argument('-j', '--jobs', type=int, default=4)
    p.add_argument('--case', choices=list(CASES), default='high', help='export: reassortment rate case')
    p.add_argument('--rep', type=int, default=0, help='export: replicate')
    p.add_argument('--out', type=Path, help='export: dataset directory')
    args = p.parse_args()
    {'simulate': simulate, 'run': run, 'evaluate': evaluate, 'export': export}[args.step](args)


if __name__ == '__main__':
    sys.exit(main())
