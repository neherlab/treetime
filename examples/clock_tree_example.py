
from Bio import Phylo
import matplotlib.pyplot as plt
from treetime import ClockTree
import numpy as np
import seaborn as sns

if __name__ == '__main__':

    plt.ion()
    base_name = 'data/H3N2_NA_allyears_NA.20'
    #root_name = 'A/Hong_Kong/JY2/1968|CY147440|1968|Hong_Kong||H3N2/8-1416'
    root_name = 'A/New_York/182/2000|CY001279|02/18/2000|USA|99_00|H3N2/1-1409'
    with open(base_name+'.metadata.csv') as date_file:
        dates = {}
        for line in date_file:
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue

    tree = Phylo.read(base_name + ".nwk", 'newick')
    tree.root_with_outgroup([n for n in tree.get_terminals() if n.name==root_name][0])
    myTree = ClockTree(gtr='Jukes-Cantor', tree = tree,
                        aln = base_name+'.fasta', verbose = 6, dates = dates)

    myTree.optimize_seq_and_branch_len(prune_short=True)
    # fix clock_rate -- to test. Do inference both for joint and marginal reconstruction
    myTree.make_time_tree(clock_rate=0.003, time_marginal=False)
    myTree.make_time_tree(clock_rate=0.003, time_marginal=True)

    # make a figure the shows the branch length interpolators for each branch in the tree
    plt.figure()
    x = np.linspace(0,0.05,100)
    leaf_count=0
    for node in myTree.tree.find_clades(order='postorder'):
        if node.up is not None:
            plt.plot(x, node.branch_length_interpolator.prob_relative(x))
        if node.is_terminal():
            leaf_count+=1
            node.ypos = leaf_count
        else:
            node.ypos = np.mean([c.ypos for c in node.clades])
    plt.yscale('log')
    plt.xlabel("branch length")
    plt.xlabel("rel. probability density")
    plt.ylim([0.01,1.2])

    # make a figure that shows the time scaled tree and the probability distributions of the
    # node positions. In addition, the branch length probabilities are added to the tree.
    fig, axs = plt.subplots(2,1, sharex=True, figsize=(8,12))
    x = np.linspace(-0.1,0.05,1000)+ myTree.tree.root.time_before_present
    Phylo.draw(tree, axes=axs[0], show_confidence=False)
    offset = myTree.tree.root.time_before_present + myTree.tree.root.branch_length
    cols = sns.color_palette()
    depth = myTree.tree.depths()
    for ni,node in enumerate(myTree.tree.find_clades()):
        if (not node.is_terminal()):
            axs[1].plot(offset-x, node.marginal_pos_LH.prob_relative(x), '-', c=cols[ni%len(cols)])
        if node.up is not None:
            x_branch = np.linspace(depth[node]-2*node.branch_length-0.005,depth[node],100)
            axs[0].plot(x_branch, node.ypos - 0.7*node.branch_length_interpolator.prob_relative(depth[node]-x_branch), '-', c=cols[ni%len(cols)])
    axs[1].set_yscale('log')
    axs[1].set_ylim([0.01,1.2])
    axs[0].set_xlabel('')
    plt.tight_layout()
