from __future__ import print_function, division
from treetime import TreeTime
import numpy as np
from scipy import optimize as sciopt
import matplotlib.pyplot as plt
from matplotlib import cm
from Bio import Phylo
try:
    import seaborn as sns
    sns.set_style('whitegrid')
except:
    pass

if __name__ == '__main__':

    # load data and parse dates
    plt.ion()
    base_name = 'data/H3N2_NA_allyears_NA.20'
    import datetime
    from treetime.utils import numeric_date
    with open(base_name+'.metadata.csv') as date_file:
        dates = {}
        for line in date_file:
            if line[0]=='#':
                continue
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue

    # instantiate treetime
    tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                          aln = base_name+'.fasta', verbose = 0, dates = dates)

    # inititally the root if the tree is a mess:
    fig, axs = plt.subplots(1,2, figsize=(12,6))
    Phylo.draw(tt.tree, show_confidence=False, axes=axs[0], label_func =lambda x:x.name.split('|')[0])
    tt.plot_root_to_tip(ax=axs[-1])

    # lets reroot: we now have a positve correlation of root-to-tip distance with sampling date
    tt.reroot(root="best")
    fig, axs = plt.subplots(1,2, figsize=(12,6))
    Phylo.draw(tt.tree, show_confidence=False, axes=axs[0], label_func =lambda x:x.name.split('|')[0])
    tt.plot_root_to_tip(ax=axs[-1])

    # rerooting can be done along with the tree time inference
    tt.run(root="best")
    # if more complicated models (relaxed clocks, coalescent models) are to be used
    # or you want to resolve polytomies, treetime needs to be run for
    # several iterations, for example as
    # tt.run(root="best", resolve_polytomies=True, max_iter=2)

    # each node is now at a position that correspond to the given or inferred date
    # the units of branch length are still clock rate.
    print("clock rate: %1.5f"%tt.date2dist.slope)
    fig, axs = plt.subplots(1,2, figsize=(12,6))
    Phylo.draw(tt.tree, label_func=lambda x:'', show_confidence=False, axes=axs[0])
    axs[1].set_title("units are substitutions")
    # we can convert the branch length to units in years and redraw
    tt.branch_length_to_years()
    Phylo.draw(tt.tree, label_func=lambda x:'', show_confidence=False, axes=axs[1])
    axs[1].set_title("units are years")

    # treetime implements a convenience function to plot timetrees
    from treetime.io import plot_vs_years
    plot_vs_years(tt, label_func=lambda x:"", show_confidence=False)


