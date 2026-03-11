from __future__ import print_function, division
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from treetime import TreeAnc
import treetime

if __name__ == '__main__':
    plt.ion()

    # path to data
    base_name = '../data/h3n2_na/h3n2_na_200'

    # instantiate treetime
    myTree = TreeAnc(gtr='Jukes-Cantor', tree=base_name+".nwk", aln=base_name+'.fasta', verbose = 0)

    # infer a GTR model and ML ancestral sequences, either marginally or jointly most likely
    myTree.infer_ancestral_sequences(infer_gtr=True, marginal=True)

    # lets examine the properties of a node in the tree after ancestral inference
    node = myTree.tree.get_nonterminals()[7]
    # each node now has an inferred sequence
    if treetime.version<"0.7":
        print("\nthe inferred sequences is an array of characters:", node.sequence)
    else:
        print("\nthe inferred sequences is an array of characters:", myTree.sequence(node, as_string=False))

    # in addition, each node of the tree now has an mutation object attached
    # note that the mutation numbering starts at 0 rather than 1
    print("\nmutations on the branch leading to node %s:"%node.name, node.mutations)

    # we can readily verify these mutations by checking the inferred sequences
    if node.mutations:
        mut = node.mutations[0]
        if treetime.version<"0.7":
            print("\nmutation %s%d%s corresponds to"%mut,
                  "parent state: %s, child state %s\n\n"%(node.up.sequence[mut[1]], node.sequence[mut[1]]))
        else:
            print("\nmutation %s%d%s corresponds to"%mut,
                  "parent state: %s, child state %s\n\n"%(myTree.sequence(node.up)[mut[1]], myTree.sequence(node)[mut[1]]))


    # plot the tree and label branches by mutation
    fig, ax = plt.subplots(1,1,figsize=(18,13))
    ax.set_title("branches annotated by inferred mutations", fontsize=18)
    ax.set_axis_off()
    Phylo.draw(myTree.tree, label_func = lambda x:"", axes=ax,
        branch_labels=lambda x:", ".join(map(lambda m: "%s%d%s"%m, x.mutations[:3]))
                               +('...' if len(x.mutations)>3 else ""))

    # finally, we print the inferred GTR model
    print(myTree.gtr)
