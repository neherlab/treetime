from __future__ import print_function, division
import numpy as np
from treetime import TreeTime
import datetime
from Bio import Phylo, AlignIO
import matplotlib.pyplot as plt
plt.ion()



if __name__=='__main__':
    base_name = 'data/H3N2_NA_allyears_NA.20'
    with open(base_name+'.metadata.csv') as date_file:
        dates = {}
        for line in date_file:
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue

    # instantiate tree time with gtr, tree, alignment, and dates
    myTree = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                        aln = base_name+'.fasta', verbose = 3, dates = dates)

    # run with rerooting to best root, with a strict clock, but
    # resolving polytomies. Coalescent timescale is set to 0.01 hamming
    # distance and the optimization is run for a maximum of two cycles.
    myTree.run(root='best', relaxed_clock=False,  resolve_polytomies=True,
               Tc=0.01, max_iter=2)

    # draw tree with branch length units
    Phylo.draw(myTree.tree, show_confidence=False, label_func = lambda x:'')

    # draw tree with branch years as units
    from treetime.io import plot_vs_years
    plot_vs_years(myTree, show_confidence=False, label_func = lambda x:'')

    r2tip = np.array([[n.numdate, n.dist2root] for n in myTree.tree.get_terminals()])
    r2int = np.array([[n.numdate, n.dist2root] for n in myTree.tree.get_nonterminals()])
    plt.figure()
    plt.scatter(r2tip[:,0], r2tip[:,1], c='g', label='terminal nodes', s=30)
    plt.scatter(r2int[:,0], r2int[:,1], c='b', label='internal nodes', s=30)
    plt.ylabel('root to node distance')
    plt.xlabel('date')
    plt.legend(loc=2)


