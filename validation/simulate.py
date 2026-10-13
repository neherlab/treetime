"""Simulate one coalescent tree and evolve sequences along it for every (L, mu) combination.

Usage: python simulate.py --coal exp --rep 3 --out data [--tips-per-sample 10]
"""

import argparse
import copy
import json
import lzma
from pathlib import Path

import numpy as np
from Bio import Phylo
from treetime.gtr import GTR
from treetime.seqgen import SeqGen

from coalescent import PRESENT, SCHEDULES, simulate


def random_gtr(mu, rng):
    """GTR with random pi and symmetric W, normalised so that mu is substitutions per site and year."""
    W = rng.random((4, 4))
    W = W + W.T
    np.fill_diagonal(W, 0)
    gtr = GTR.custom(mu=1.0, pi=rng.random(4), W=W, alphabet='nuc_nogap')
    gtr.mu = mu  # assign_rates multiplies mu by the raw average rate; W is normalised now
    return gtr


def write_fasta_xz(path, nodes):
    with lzma.open(path, 'wt') as fh:
        for n in nodes:
            fh.write(f'>{n.name}\n{"".join(n.ancestral_sequence)}\n')


def write_dates(path, nodes, tau):
    with open(path, 'w') as fh:
        fh.write('name\tdate\ttau\n')
        for n in nodes:
            fh.write(f'{n.name}\t{PRESENT - tau[n]:.8f}\t{tau[n]:.8f}\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--coal', choices=SCHEDULES, required=True)
    p.add_argument('--rep', type=int, required=True)
    p.add_argument('--out', type=Path, default=Path('data'))
    p.add_argument('--lengths', type=int, nargs='+', default=[1000, 10000])
    p.add_argument('--mus', type=float, nargs='+', default=[0.01, 0.003, 0.001, 0.0003])
    p.add_argument('--seed', type=int, default=2025)
    p.add_argument('--tips-per-sample', type=int, default=10, help='tips sampled at each of the 40 time points')
    args = p.parse_args()

    key = [args.seed, list(SCHEDULES).index(args.coal), args.rep]
    if args.tips_per_sample != 10:  # the 10-tip data keeps its original seeds
        key.append(args.tips_per_sample)
    outdir = args.out / args.coal / f'rep{args.rep:02d}'
    outdir.mkdir(parents=True, exist_ok=True)

    tree, tau = simulate(SCHEDULES[args.coal](), np.random.default_rng(key), n_per_batch=args.tips_per_sample)
    tree.ladderize()
    Phylo.write(tree, outdir / 'tree.nwk', 'newick', format_branch_length='%1.8f')
    write_dates(outdir / 'tip_dates.tsv', tree.get_terminals(), tau)
    write_dates(outdir / 'node_dates.tsv', tree.get_nonterminals(), tau)

    for L in args.lengths:
        for mu in args.mus:
            sub = outdir / f'L{L}_mu{mu:g}'
            if (sub / 'gtr.json').exists():  # written last, marks a complete set
                continue
            rng = np.random.default_rng([*key, L, round(mu * 1e6)])
            gtr = random_gtr(mu, rng)
            sim = SeqGen(L, tree=copy.deepcopy(tree), gtr=gtr, rng_seed=int(rng.integers(2**32)), verbose=0)
            sim.evolve()
            sub.mkdir(exist_ok=True)
            write_fasta_xz(sub / 'tips.fasta.xz', sim.tree.get_terminals())
            write_fasta_xz(sub / 'internal.fasta.xz', sim.tree.get_nonterminals())
            with open(sub / 'gtr.json', 'w') as fh:
                json.dump(
                    {'alphabet': 'ACGT', 'mu': gtr.mu, 'pi': gtr.Pi.tolist(), 'W': gtr.W.tolist(), 'L': L}, fh, indent=1
                )
            print(f'{outdir} L={L} mu={mu:g} done', flush=True)


if __name__ == '__main__':
    main()
