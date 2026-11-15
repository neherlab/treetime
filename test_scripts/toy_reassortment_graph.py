"""Generate the toy reassortment graph in data/toy/reassortment/.

A hand-made graph with three segments. Every edge carries a set of segments
(its scope). For each segment, the edges that carry it form one tree (the
local tree of that segment). The script writes, per segment:

- the segment tree as Newick, with every node of the segment's local tree,
  including nodes that the segment only passes through (one child)
- an alignment of the leaves, simulated along the local tree under
  Jukes-Cantor at rate 1 (branch lengths are substitutions per site)
- the true sequences of all nodes of the local tree

and for the whole graph an edge table with the scopes, node times, and leaf
dates.

Usage: python3 test_scripts/toy_reassortment_graph.py [output_dir]
"""

import sys
from pathlib import Path

import numpy as np

SEED = 20261115
SEGMENTS = {"seg1": 600, "seg2": 500, "seg3": 400}
ALPHABET = np.array(list("ACGT"))

# Node times in substitutions per site from the oldest root.
TIMES = {
    "R": 0.000,
    "R3": 0.005,
    "N": 0.015,
    "X": 0.020,
    "Y": 0.030,
    "K": 0.040,
    "M": 0.045,
    "P": 0.050,
    "U": 0.055,
    "Q": 0.060,
    "Z": 0.065,
    "V": 0.070,
    "A": 0.100,
    "B": 0.095,
    "C": 0.090,
    "D": 0.100,
    "E": 0.085,
    "F": 0.095,
    "G": 0.100,
    "H": 0.090,
    "I": 0.095,
    "J": 0.100,
}

ALL = ("seg1", "seg2", "seg3")

# (parent, child, segments carried). Child order follows this list.
EDGES = [
    ("R", "X", ("seg1", "seg2")),
    ("R", "Y", ("seg1", "seg2")),
    ("R3", "N", ("seg3",)),
    ("R3", "X", ("seg3",)),
    ("N", "Y", ("seg3",)),
    ("N", "V", ("seg3",)),
    ("X", "P", ("seg1", "seg3")),
    ("X", "M", ALL),
    ("Y", "K", ALL),
    ("Y", "V", ("seg1",)),
    ("Y", "Z", ALL),
    ("K", "U", ALL),
    ("K", "P", ("seg2",)),
    ("M", "Q", ALL),
    ("M", "V", ("seg2",)),
    ("P", "A", ALL),
    ("P", "B", ALL),
    ("Q", "C", ALL),
    ("Q", "D", ALL),
    ("U", "E", ALL),
    ("U", "F", ALL),
    ("V", "G", ALL),
    ("V", "H", ALL),
    ("Z", "I", ALL),
    ("Z", "J", ALL),
]


def year(t):
    return 2000.0 + 200.0 * t


def local_tree(segment):
    children = {}
    parents = {}
    for parent, child, scope in EDGES:
        if segment in scope:
            if child in parents:
                raise ValueError(f"{child} has two parents in {segment}")
            parents[child] = parent
            children.setdefault(parent, []).append(child)
    nodes = set(children) | set(parents)
    roots = sorted(n for n in nodes if n not in parents)
    return roots, children


def check_graph():
    for parent, child, _ in EDGES:
        if TIMES[child] <= TIMES[parent]:
            raise ValueError(f"edge {parent}->{child} goes back in time")
    leaves = {c for _, c, _ in EDGES} - {p for p, _, _ in EDGES}
    for segment in SEGMENTS:
        roots, children = local_tree(segment)
        if len(roots) != 1:
            raise ValueError(f"{segment} has roots {roots}")
        reached = set()
        stack = list(roots)
        while stack:
            node = stack.pop()
            reached.add(node)
            stack.extend(children.get(node, []))
        if not leaves <= reached:
            raise ValueError(f"{segment} misses leaves {sorted(leaves - reached)}")


def newick(node, children, parent_time=None):
    kids = children.get(node, [])
    inner = "(" + ",".join(newick(c, children, TIMES[node]) for c in kids) + ")" if kids else ""
    length = "" if parent_time is None else f":{TIMES[node] - parent_time:.6f}"
    return f"{inner}{node}{length}"


def evolve(rng, seq, t):
    p_change = 0.75 * (1.0 - np.exp(-4.0 / 3.0 * t))
    change = rng.random(seq.size) < p_change
    shift = rng.integers(1, 4, size=seq.size)
    return np.where(change, (seq + shift) % 4, seq)


def simulate(rng, segment, length):
    roots, children = local_tree(segment)
    seqs = {}
    for root in roots:
        seqs[root] = rng.integers(0, 4, size=length)
        stack = [root]
        while stack:
            node = stack.pop()
            for child in children.get(node, []):
                seqs[child] = evolve(rng, seqs[node], TIMES[child] - TIMES[node])
                stack.append(child)
    return roots, children, seqs


def write_fasta(path, seqs, names):
    with open(path, "w") as f:
        for name in names:
            f.write(f">{name}\n{''.join(ALPHABET[seqs[name]])}\n")


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/toy/reassortment")
    out.mkdir(parents=True, exist_ok=True)
    check_graph()
    rng = np.random.default_rng(SEED)
    leaves = sorted({c for _, c, _ in EDGES} - {p for p, _, _ in EDGES})
    for segment, length in SEGMENTS.items():
        roots, children, seqs = simulate(rng, segment, length)
        (out / f"tree_{segment}.nwk").write_text(newick(roots[0], children) + ";\n")
        write_fasta(out / f"aln_{segment}.fasta", seqs, leaves)
        internal = sorted((n for n in seqs if n not in leaves), key=lambda n: (TIMES[n], n))
        write_fasta(out / f"truth_{segment}.fasta", seqs, internal + leaves)

    with open(out / "graph.tsv", "w") as f:
        f.write("parent\tchild\tlength\tsegments\n")
        for parent, child, scope in EDGES:
            f.write(f"{parent}\t{child}\t{TIMES[child] - TIMES[parent]:.6f}\t{','.join(scope)}\n")

    with open(out / "nodes.tsv", "w") as f:
        f.write("name\ttime\tdate\n")
        for name, t in sorted(TIMES.items(), key=lambda kv: (kv[1], kv[0])):
            f.write(f"{name}\t{t:.6f}\t{year(t):.2f}\n")

    with open(out / "metadata.tsv", "w") as f:
        f.write("name\tdate\n")
        for name in leaves:
            f.write(f"{name}\t{year(TIMES[name]):.2f}\n")


if __name__ == "__main__":
    main()
