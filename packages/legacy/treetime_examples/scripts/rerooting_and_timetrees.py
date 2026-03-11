from __future__ import print_function, division
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib import cm
try:
    import seaborn as sns
    sns.set_style('whitegrid')
except:
    print ("Seaborn not found. Default style will be used for the plots")

from treetime import TreeTime
from treetime.utils import parse_dates

def format_axes(fig, axs):
    axs[0].set_axis_off()
    axs[1].tick_params(labelsize=14)
    axs[1].set_ylabel("root-to-tip distance", fontsize=16)
    axs[1].set_xlabel("date", fontsize=16)
    fig.tight_layout()


if __name__ == '__main__':

    # load data and parse dates
    plt.ion()
    base_name = '../data/h3n2_na/h3n2_na_20'

    dates = parse_dates(base_name+'.metadata.csv')
    tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                  aln = base_name+'.fasta', verbose = 1, dates = dates)

    # inititally the root if the tree is a mess:
    fig, axs = plt.subplots(1,2, figsize=(18,9))
    axs[0].set_title("Arbitrarily rooted tree", fontsize=18)
    axs[1].set_title("Inverse divergence-time relationship", fontsize=18)
    Phylo.draw(tt.tree, show_confidence=False, axes=axs[0], label_func=lambda x:x.name.split('|')[0] if x.is_terminal() else "")
    tt.plot_root_to_tip(ax=axs[-1])
    format_axes(fig, axs)

    # lets reroot: we now have a positve correlation of root-to-tip distance with sampling date
    tt.reroot(root="least-squares")
    fig, axs = plt.subplots(1,2, figsize=(18,9))
    axs[0].set_title("Tree rerooted by treetime", fontsize=18)
    axs[1].set_title("Optimal divergence-time relationship", fontsize=18)
    Phylo.draw(tt.tree, show_confidence=False, axes=axs[0],
               label_func=lambda x:x.name.split('|')[0] if x.is_terminal() else "")
    tt.plot_root_to_tip(ax=axs[-1])
    format_axes(fig, axs)

    # rerooting can be done along with the tree time inference
    tt.run(root="best", branch_length_mode='input', max_iter=2)
    # if more complicated models (relaxed clocks, coalescent models) are to be used
    # or you want to resolve polytomies, treetime needs to be run for
    # several iterations, for example as
    # tt.run(root="best", resolve_polytomies=True, max_iter=2)

    # each node is now at a position that correspond to the given or inferred date
    # the units of branch length are still clock rate.
    print("clock rate: %1.5f"%tt.date2dist.clock_rate)
    fig, axs = plt.subplots(1,2, figsize=(18,9))
    Phylo.draw(tt.tree, label_func=lambda x:'', show_confidence=False, axes=axs[0])
    axs[0].set_title("Tree: units are substitutions", fontsize=18)

    # we can convert the branch length to units in years and redraw
    tt.branch_length_to_years()
    Phylo.draw(tt.tree, label_func=lambda x:'', show_confidence=False, axes=axs[1])
    axs[1].set_title("Tree: units are years", fontsize=18)
    axs[0].tick_params(labelsize=14)
    axs[1].tick_params(labelsize=14)
    fig.tight_layout()

    # treetime implements a convenience function to plot timetrees
    from treetime.treetime import plot_vs_years
    plot_vs_years(tt, label_func=lambda x:"", show_confidence=False)

