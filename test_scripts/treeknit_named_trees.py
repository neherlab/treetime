"""Convert TreeKnit output into segment trees whose shared internal nodes have equal names.

TreeKnit writes one resolved tree per segment, with internal names local to each tree, and for every pair of
segments the maximum compatible clades (MCCs): sets of leaves whose subtree topology both trees share. This script

1. assigns each internal node of a tree, for each partner segment, the MCC it belongs to (the rule of v0
   `treetime/arg.py:assign_mccs`: Fitch on the leaf MCCs, then inheritance from the parent where possible)
2. matches, for each segment pair and each MCC, the node of one tree with the node of the other tree that is the
   most recent common ancestor (MRCA) of the same MCC leaves below it
3. joins the matches of all pairs into node identities (union-find) and reports identities that would hold two
   nodes of one tree (conflicts between pairwise MCCs); such identities are split back into separate nodes
4. splits identities whose segments order them in opposite ways (a cycle in the merged graph) back into
   separate nodes, one cycle at a time
5. writes the trees with one name per identity

Usage: python3 -I test_scripts/treeknit_named_trees.py <treeknit_dir> <output_dir>
"""

import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import Phylo

SEGMENTS = ["ha", "na", "pb1", "pb2"]


def read_tree(directory, segment):
    tree = Phylo.read(directory / f"tree_{segment}_resolved.nwk", "newick")
    parents = {}
    for clade in tree.find_clades(order="preorder"):
        for child in clade.clades:
            parents[child] = clade
    leaves = {}
    for clade in tree.find_clades(order="postorder"):
        leaves[clade] = frozenset([clade.name]) if clade.is_terminal() else frozenset().union(*(leaves[c] for c in clade.clades))
    return tree, parents, leaves


def assign_mccs(tree, parents, leaf_mcc):
    child_mccs = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            child_mccs[clade] = {leaf_mcc[clade.name]}
            continue
        sets = [child_mccs[c] for c in clade.clades]
        common = set.intersection(*sets)
        child_mccs[clade] = common if common else set.union(*sets)
    mcc = {}
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            mcc[clade] = leaf_mcc[clade.name]
        elif clade not in parents:
            common = set.intersection(*(child_mccs[c] for c in clade.clades))
            mcc[clade] = min(common) if common else None
        elif mcc[parents[clade]] in child_mccs[clade]:
            mcc[clade] = mcc[parents[clade]]
        elif len(child_mccs[clade]) == 1:
            mcc[clade] = next(iter(child_mccs[clade]))
        else:
            mcc[clade] = None
    return mcc


def mrca(tree_data, leaf_set):
    tree, parents, leaves = tree_data
    node = tree.find_any(name=next(iter(leaf_set)))
    while not leaf_set <= leaves[node]:
        node = parents[node]
    return node


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        self.parent[self.find(a)] = self.find(b)


def find_cycle(trees, names):
    children = defaultdict(set)
    for segment, (tree, parents, _) in trees.items():
        for child, parent in parents.items():
            children[node_name(names, segment, parent)].add(node_name(names, segment, child))
    state = {}
    for start in list(children):
        if start in state:
            continue
        state[start] = "open"
        path, stack = [start], [iter(children[start])]
        while stack:
            child = next(stack[-1], None)
            if child is None:
                state[path.pop()] = "done"
                stack.pop()
            elif state.get(child) == "open":
                return path[path.index(child):]
            elif child not in state:
                state[child] = "open"
                path.append(child)
                stack.append(iter(children[child]))
    return None


def node_name(names, segment, clade):
    return names.get((segment, clade), clade.name)


def break_cycles(trees, names):
    splits = []
    while cycle := find_cycle(trees, names):
        shared = next(name for name in cycle if name.startswith("N"))
        for (segment, clade), name in list(names.items()):
            if name == shared:
                names[(segment, clade)] = f"{segment}_{clade.name}"
        splits.append(shared)
    return splits


def main():
    directory, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    trees = {segment: read_tree(directory, segment) for segment in SEGMENTS}
    pairs = json.load(open(directory / "MCCs.json"))["MCC_dict"]

    identities = UnionFind()
    stats = {}
    for entry in pairs.values():
        a, b = (label.removeprefix("tree_") for label in entry["trees"])
        leaf_mcc = {leaf: index for index, mcc in enumerate(entry["mccs"]) for leaf in mcc}
        mccs = {segment: assign_mccs(trees[segment][0], trees[segment][1], leaf_mcc) for segment in (a, b)}
        matched = unmatched = 0
        for node_b, m in mccs[b].items():
            if m is None or node_b.is_terminal():
                continue
            restricted = trees[b][2][node_b] & frozenset(entry["mccs"][m])
            node_a = mrca(trees[a], restricted)
            if mccs[a].get(node_a) == m and not node_a.is_terminal() and mrca(trees[b], restricted) is node_b:
                identities.union((a, node_a), (b, node_b))
                matched += 1
            else:
                unmatched += 1
        stats[f"{a}-{b}"] = {"mccs": len(entry["mccs"]), "matched_nodes": matched, "unmatched_nodes": unmatched}

    classes = defaultdict(list)
    for segment, (tree, _, _) in trees.items():
        for clade in tree.get_nonterminals():
            classes[identities.find((segment, clade))].append((segment, clade))
    conflicts = [members for members in classes.values() if len({s for s, _ in members}) < len(members)]
    sizes = Counter(len({s for s, _ in members}) for members in classes.values() if members not in conflicts)

    names = {}
    for index, members in enumerate(sorted(classes.values(), key=lambda m: (m[0][0], m[0][1].name))):
        if members in conflicts:
            for segment, clade in members:
                names[(segment, clade)] = f"{segment}_{clade.name}"
        else:
            for segment, clade in members:
                names[(segment, clade)] = f"N{index:05d}"
    cycle_splits = break_cycles(trees, names)
    for segment, (tree, _, _) in trees.items():
        for clade in tree.get_nonterminals():
            clade.name = names[(segment, clade)]
            clade.confidence = None
        Phylo.write(tree, out / f"tree_{segment}.nwk", "newick", format_branch_length="%.8f")

    report = {
        "pairs": stats,
        "identities_by_number_of_segments": dict(sorted(sizes.items())),
        "conflicting_identities": len(conflicts),
        "nodes_in_conflicts": sum(len(m) for m in conflicts),
        "identities_split_to_break_cycles": cycle_splits,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
