"""Divergence trees (substitutions per site) derived from the true time tree of a simulated sequence set.

exact      branch length mu * t
realized   true parent-child differences / L
rate       realized * lognormal(0, 0.25) per branch
"""

import json
import lzma
import zlib

import numpy as np
from Bio import Phylo, SeqIO

RATE_SIGMA = 0.25


def read_fasta_xz(path):
    with lzma.open(path, 'rt') as fh:
        return {r.id: np.frombuffer(str(r.seq).encode(), 'S1') for r in SeqIO.parse(fh, 'fasta')}


def divergence_tree(set_dir, kind):
    """Write (once) and return the path of the divergence tree `kind` for the sequence set in `set_dir`."""
    path = set_dir / 'inputs' / f'{kind}.nwk'
    if path.exists():
        return path
    rng = np.random.default_rng(zlib.crc32(f'{"/".join(set_dir.parts[-3:])}:{kind}'.encode()))
    tree = Phylo.read(set_dir.parent / 'tree.nwk', 'newick')
    gtr = json.loads((set_dir / 'gtr.json').read_text())

    if kind == 'exact':
        for c in tree.find_clades():
            c.branch_length = (c.branch_length or 0.0) * gtr['mu']
    else:
        seqs = read_fasta_xz(set_dir / 'tips.fasta.xz') | read_fasta_xz(set_dir / 'internal.fasta.xz')
        for p in tree.get_nonterminals():
            for c in p:
                c.branch_length = np.count_nonzero(seqs[c.name] != seqs[p.name]) / gtr['L']
                if kind == 'rate':
                    c.branch_length *= rng.lognormal(0.0, RATE_SIGMA)
        tree.root.branch_length = 0.0

    path.parent.mkdir(exist_ok=True)
    tmp = path.with_suffix('.tmp')
    Phylo.write(tree, tmp, 'newick', format_branch_length='%1.10g')
    tmp.rename(path)
    return path
