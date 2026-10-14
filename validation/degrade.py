"""Make a messy copy of one simulated replicate: phylogenetic gaps, ambiguous bases, runs of N, imprecise tip dates.

Usage: python degrade.py --coal exp --rep 3 [--data data] [--out data_messy]

The copy has the layout of the input (see README), with truth and observations kept apart:
  tree.nwk, node_dates.tsv, tip_dates.tsv     unchanged truth
  tip_dates_observed.tsv                      name, date: some dates missing (empty), some only to the month (YYYY-MM-XX)
  L<L>_mu<mu>/
    tips.fasta.xz, internal.fasta.xz          truth with gaps: a gap on a branch is inherited by every node below it
    tips_observed.fasta.xz                    tips with gaps, ambiguous bases (IUPAC codes containing the true base)
                                              and runs of N
    gaps.tsv                                  node, start of each 3-base gap
    gtr.json                                  unchanged, written last
run.py hands the *_observed files to TreeTime when they exist.
"""

import argparse
import datetime
import json
import lzma
import shutil
import zlib
from pathlib import Path

import numpy as np
from Bio import Phylo

from inputs import read_fasta_xz

GAP_LENGTH = 3
GAPS_PER_KB = 3  # gap events per 1000 sites, on branches drawn in proportion to their length in time
AMBIGUOUS_RATE = 0.002  # per tip and non-gap site
N_RUN_TIP_FRACTION = 0.05  # tips with one run of N
N_RUN_LENGTH = (5, 100)
MISSING_DATE_FRACTION = 0.03
MONTH_DATE_FRACTION = 0.05
IUPAC = {'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG'}
CODES_WITH = {b: np.array([c for c, s in IUPAC.items() if b in s], 'S1') for b in 'ACGT'}


def rng_for(*key):
    return np.random.default_rng(zlib.crc32(':'.join(map(str, key)).encode()))


def month_only(numdate):
    year = int(numdate)
    start = datetime.date(year, 1, 1)
    days = (datetime.date(year + 1, 1, 1) - start).days
    d = start + datetime.timedelta(days=(numdate - year) * days)
    return f'{d.year}-{d.month:02d}-XX'


def degrade_dates(src, dst, rng):
    with open(src) as fh, open(dst, 'w') as out:
        next(fh)
        out.write('name\tdate\n')
        for line in fh:
            name, date, _ = line.rstrip('\n').split('\t')
            u = rng.random()
            if u < MISSING_DATE_FRACTION:
                date = ''
            elif u < MISSING_DATE_FRACTION + MONTH_DATE_FRACTION:
                date = month_only(float(date))
            out.write(f'{name}\t{date}\n')


def add_gaps(tree, seqs, L, rng):
    """Delete GAP_LENGTH sites on random branches, in the node and all its descendants. Returns (node, start) pairs."""
    branches = [c for c in tree.find_clades() if c is not tree.root]
    p = np.array([c.branch_length for c in branches])
    events = []
    for i in rng.choice(len(branches), size=max(1, round(GAPS_PER_KB * L / 1000)), p=p / p.sum()):
        start = int(rng.integers(L - GAP_LENGTH + 1))
        for c in branches[i].find_clades():
            seqs[c.name][start : start + GAP_LENGTH] = b'-'
        events.append((branches[i].name, start))
    return events


def observe(seq, rng):
    """Tip sequence as sequenced: ambiguous codes containing the true base, and sometimes a run of N."""
    seq = seq.copy()
    for i in np.flatnonzero((rng.random(len(seq)) < AMBIGUOUS_RATE) & (seq != b'-')):
        seq[i] = rng.choice(CODES_WITH[seq[i].decode()])
    if rng.random() < N_RUN_TIP_FRACTION:
        n = int(rng.integers(N_RUN_LENGTH[0], N_RUN_LENGTH[1] + 1))
        start = int(rng.integers(len(seq) - n + 1))
        seq[start : start + n] = b'N'
    return seq


def write_fasta_xz(path, seqs, names):
    with lzma.open(path, 'wt') as fh:
        for n in names:
            fh.write(f'>{n}\n{seqs[n].tobytes().decode()}\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--coal', required=True)
    p.add_argument('--rep', type=int, required=True)
    p.add_argument('--data', type=Path, default=Path('data'))
    p.add_argument('--out', type=Path, help='default: <data>_messy')
    args = p.parse_args()

    rep = f'rep{args.rep:02d}'
    src = args.data / args.coal / rep
    dst = (args.out or args.data.with_name(args.data.name + '_messy')) / args.coal / rep
    dst.mkdir(parents=True, exist_ok=True)
    for f in ['tree.nwk', 'node_dates.tsv', 'tip_dates.tsv']:
        shutil.copy(src / f, dst / f)
    degrade_dates(src / 'tip_dates.tsv', dst / 'tip_dates_observed.tsv', rng_for(args.data.name, args.coal, rep))
    tree = Phylo.read(src / 'tree.nwk', 'newick')
    tips = [c.name for c in tree.get_terminals()]
    internal = [c.name for c in tree.get_nonterminals()]

    for set_src in sorted(d for d in src.glob('L*') if (d / 'gtr.json').exists()):
        set_dst = dst / set_src.name
        if (set_dst / 'gtr.json').exists():  # written last, marks a complete set
            continue
        set_dst.mkdir(exist_ok=True)
        rng = rng_for(args.data.name, args.coal, rep, set_src.name)
        seqs = {k: v.copy() for k, v in (read_fasta_xz(set_src / 'tips.fasta.xz') | read_fasta_xz(set_src / 'internal.fasta.xz')).items()}
        L = len(seqs[tips[0]])
        gaps = add_gaps(tree, seqs, L, rng)
        (set_dst / 'gaps.tsv').write_text('node\tstart\n' + ''.join(f'{n}\t{s}\n' for n, s in gaps))
        write_fasta_xz(set_dst / 'tips.fasta.xz', seqs, tips)
        write_fasta_xz(set_dst / 'internal.fasta.xz', seqs, internal)
        write_fasta_xz(set_dst / 'tips_observed.fasta.xz', {n: observe(seqs[n], rng) for n in tips}, tips)
        shutil.copy(set_src / 'gtr.json', set_dst / 'gtr.json')
        print(f'{set_dst} done', flush=True)


if __name__ == '__main__':
    main()
