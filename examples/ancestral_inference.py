from __future__ import print_function, division
from treetime import TreeAnc
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from Bio import Phylo

if __name__ == '__main__':
    plt.ion()

    # load data
    base_name = 'data/H3N2_NA_allyears_NA.20'
    T = Phylo.read(base_name+".nwk", "newick")
    T.root_with_outgroup(T.find_clades(lambda x:x.name.startswith('A/Scot')).next())

    # instantiate treetime
    myTree = TreeAnc(gtr='Jukes-Cantor', tree = T, aln = base_name+'.fasta', verbose = 0)

    # infer optimize branch length, infer a GTR model, and reoptimize branch length
    myTree.optimize_sequences_and_branch_length(infer_gtr=True)

    # lets examine the properties of a node in the tree after ancestral inference
    node = myTree.tree.get_nonterminals()[7]
    # each node now has an inferred sequence
    print("\nthe inferred sequences is an array of states:", node.sequence)

    # in addition, each node of the tree now has an mutation object attached
    # note that the mutation numbering starts at 0 rather than 1
    print("\nmutations of node %s:"%node.name, node.mutations)

    # we can readily verify these mutations by checking the inferred sequences
    mut = node.mutations[0]
    print("\nmutation %s%d%s corresponds to"%mut,
          "parent state: %s, child state %s\n\n"%(node.up.sequence[mut[1]], node.sequence[mut[1]]))


    # plot the tree and label branches by mutation
    plt.figure(figsize=(18,13))
    plt.title("branches annotated by inferred mutations", fontsize=18)
    ax=plt.subplot(111)
    ax.set_axis_off()
    plt.title
    Phylo.draw(myTree.tree, label_func = lambda x:"", axes=ax,
        branch_labels=lambda x:", ".join(map(lambda m: "%s%d%s"%m, x.mutations[:3]))
                               +('...' if len(x.mutations)>3 else ""))

    # finally, we print the inferred GTR model
    print(myTree.gtr)
