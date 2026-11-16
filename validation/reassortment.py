"""Coalescent with reassortment: the ancestral reassortment graph (ARG) of a segmented genome.

Simulated backwards in time, as in the CoalRe model (Mueller, Kistler and Bedford 2022, PLoS Comp Biol 18:e1009985):
each lineage carries the set of segments whose ancestry it holds. A pair of lineages merges at rate 1/Tc(tau), and
the merged lineage carries the union of their segments. A lineage that carries more than one segment reassorts at
rate `rho`: each of its segments goes to one of two parent lineages with probability 1/2; a split that sends every
segment to the same parent leaves the lineage unchanged. A segment stops being tracked once one lineage carries all
its ancestry, so every segment has one root.

The ARG is returned as nodes with times and as edges (parent, child, segments). The local tree of a segment is the
set of edges that carry it; it contains nodes with one child where the segment only passes through: reassortment
nodes, and merges in which only the other parent carries the segment.
"""

from dataclasses import dataclass, field

import numpy as np

from coalescent import PRESENT


@dataclass
class Arg:
    tau: dict = field(default_factory=dict)
    edges: list = field(default_factory=list)
    tips: list = field(default_factory=list)
    n_reassortments: int = 0


def simulate_arg(schedule, rng, segments, rho, n_per_batch=2, n_batches=40, interval=0.1):
    arg = Arg()
    sample_times = [i * interval for i in range(n_batches)]
    active = []
    remaining = {s: n_per_batch * n_batches for s in segments}
    carriers = dict.fromkeys(segments, 0)
    counters = {'tip': 0, 'node': 0}
    t, next_batch = 0.0, 0

    def new_node(kind, time):
        counters[kind] += 1
        name = f'{kind}_{counters[kind]:04d}'
        arg.tau[name] = time
        return name

    def retire(lineage_segments):
        done = {s for s in lineage_segments if remaining[s] == 0 and carriers[s] == 1}
        for s in done:
            carriers[s] -= 1
        return lineage_segments - done

    while True:
        t_sample = sample_times[next_batch] if next_batch < n_batches else np.inf
        k = len(active)
        coal_rate = k * (k - 1) / 2
        reass_rate = rho * sum(len(segs) > 1 for _, segs in active)
        t_coal = schedule.inv(schedule.cum(t) + rng.exponential() / coal_rate) if k > 1 else np.inf
        t_reass = t + rng.exponential() / reass_rate if reass_rate > 0 else np.inf
        t_next = min(t_sample, t_coal, t_reass)
        if t_next == np.inf:
            break
        t = t_next
        if t_next == t_sample:
            for _ in range(n_per_batch):
                tip = new_node('tip', t)
                arg.tips.append(tip)
                for s in segments:
                    remaining[s] -= 1
                    carriers[s] += 1
                active.append((tip, frozenset(segments)))
            next_batch += 1
        elif t_next == t_coal:
            i, j = sorted(rng.choice(k, size=2, replace=False), reverse=True)
            (a, sa), (b, sb) = active.pop(i), active.pop(j)
            parent = new_node('node', t)
            arg.edges += [(parent, a, sa), (parent, b, sb)]
            for s in sa & sb:
                carriers[s] -= 1
            merged = retire(sa | sb)
            if merged:
                active.append((parent, frozenset(merged)))
        else:
            weights = np.array([len(segs) > 1 for _, segs in active], dtype=float)
            i = rng.choice(len(active), p=weights / weights.sum())
            child, segs = active[i]
            side = rng.random(len(segs)) < 0.5
            if side.all() or not side.any():
                continue
            ordered = sorted(segs)
            left = frozenset(s for s, x in zip(ordered, side) if x)
            right = segs - left
            node = new_node('node', t)
            arg.n_reassortments += 1
            arg.edges.append((node, child, segs))
            active[i] = (node, left)
            active.append((node, right))
        if not active and next_batch == n_batches:
            break
    return arg


def local_tree(arg, segment):
    """(root, children, parent) of the local tree of `segment`: every node it passes through, in edge order."""
    children, parent = {}, {}
    for p, c, segs in arg.edges:
        if segment in segs:
            children.setdefault(p, []).append(c)
            parent[c] = p
    roots = [n for n in children if n not in parent]
    if len(roots) != 1:
        raise ValueError(f'segment {segment} has roots {roots}')
    return roots[0], children, parent


def collapsed(root, children, removed=None):
    """The local tree without its nodes of one child (only those in `removed`, if given)."""
    kept_children, kept_parent = {}, {}

    def below(node):
        while len(children.get(node, [])) == 1 and (removed is None or node in removed):
            node = children[node][0]
        return node

    root = below(root)
    stack = [root]
    while stack:
        node = stack.pop()
        for c in children.get(node, []):
            c = below(c)
            kept_children.setdefault(node, []).append(c)
            kept_parent[c] = node
            stack.append(c)
    return root, kept_children, kept_parent


def hidden_nodes(local_trees):
    """Nodes with one child in every local tree that holds them: reassortment nodes, and mergers of lineages that
    carry disjoint sets of segments. No segment branches there, so no tree inference recovers them."""
    unary = {}
    for _, children, parent in local_trees:
        for node in parent.keys() | children.keys():
            unary[node] = unary.get(node, True) and len(children.get(node, [])) == 1
    return {node for node, one_child in unary.items() if one_child}


def date(tau):
    return PRESENT - tau
