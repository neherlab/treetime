'''
Simple script that illustrates the usage of TreeTime. TreeTime is a subclass
of ClockTree that adds several features such as polytomy resolution, relaxed clock models,
coalescent models, rerooting, and filtering based on clock deviations.

In addition, TreeTime implements a run function that provides a unified interface to the
the different features and runs those analyses in a sensible order.
'''
from Bio import Phylo
import datetime
from treetime.utils import numeric_date, tree_layout
from treetime import TreeTime
from rerooting_and_timetrees import read_dates
import numpy as np

import matplotlib.pyplot as plt
try:
    import seaborn as sns
    sns.set_style('whitegrid')
except:
    print ("Seaborn not found. Default style will be used for the plots")
plt.ion()

if __name__ == '__main__':

    base_name = 'data/H3N2_NA_allyears_NA.20'
    dates = read_dates(base_name)

    # instantiate treetime
    myTree = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', aln = base_name+'.fasta',
                      verbose = 4, dates = dates, debug=False)

    # RUN: this example uses a fixed clock rate of 0.003, resolves polytomes, and estimates
    # confidences via a final marginal reconstruction
    myTree.run(root='clock_filter', relaxed_clock=False, max_iter=2, fixed_clock_rate=0.003,
               resolve_polytomies=True, Tc=None, n_iqd=2, time_marginal=True)



    #################
    # PLOTTING
    #################

    # draw phylogenetic tree in one panel, marginal distributions in the other
    tree_layout(myTree.tree)
    fig, axs = plt.subplots(2,1, sharex=True, figsize=(8,12))
    Phylo.draw(myTree.tree, axes=axs[0], show_confidence=False, label_func = lambda x:'')
    offset = myTree.tree.root.time_before_present + myTree.tree.root.branch_length
    cols = sns.color_palette()
    depth = myTree.tree.depths()
    x = np.linspace(-0.01, .2,1000)
    for ni,node in enumerate(myTree.tree.find_clades(order="postorder")):
        axs[1].plot(offset-x, node.marginal_pos_LH.prob_relative(x), '-', c=cols[ni%len(cols)])
        if node.up is not None:
            # add branch length distributions to tree
            x_branch = np.linspace(depth[node]-2*node.branch_length-0.005,depth[node],100)
            axs[0].plot(x_branch, node.ypos - 0.7*node.branch_length_interpolator.prob_relative(depth[node]-x_branch),
                        '-', c=cols[ni%len(cols)])
    axs[1].set_yscale('log')
    axs[1].set_ylim([0.05,1.2])
    axs[0].set_xlabel('')
    plt.tight_layout()

    # make root-to-tip plot
    myTree.plot_root_to_tip(add_internal=True, s=30)
                            # include internal nodes, set symbol size to 30