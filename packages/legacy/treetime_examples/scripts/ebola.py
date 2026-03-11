'''
this script illustrates the use of treetime to analyze viral sequences.
As an example, we use Ebola virus sequences from the 2014-2015 outbreak in
West Africa. This example contains more than 300 sequences, it will take
a few minutes to run.
'''

from __future__ import print_function, division
import sys
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
import treetime.config as ttconf

if __name__ == '__main__':


    plt.ion()
    base_name = '../data/ebola/ebola'
    dates = parse_dates(base_name+'.metadata.csv')
    # instantiate treetime
    ebola = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', precision=1,
                     aln = base_name+'.fasta', verbose = 2, dates = dates)

    # infer an ebola time tree while rerooting and resolving polytomies
    res = ebola.run(root='best', infer_gtr=True, relaxed_clock=False, max_iter=2,
                    branch_length_mode='input', n_iqd=3, resolve_polytomies=True,
                    Tc='skyline', time_marginal="assign", vary_rate=True)

    if res==ttconf.ERROR:
        sys.exit()

    # scatter root to tip divergence vs sampling date
    ebola.plot_root_to_tip(add_internal=True)
    plt.legend(loc=2)

    # rescale branch length to years and plot in axis 0
    from treetime.treetime import plot_vs_years
    fig, axs = plt.subplots(1,2, sharey=True, figsize=(12,8))
    plot_vs_years(ebola, step=1, ax=axs[1], confidence=(0.05,0.95), label_func = lambda x:"")
    axs[1].set_xlim(0, 2.5)
    axs[1].set_title("time tree")

    # assign mutation length to branch length and plot the mutation tree in axs[1]
    axs[0].set_title("mutation tree")
    for n in ebola.tree.find_clades():
        n.branch_length=n.mutation_length
    Phylo.draw(ebola.tree, label_func=lambda x:"", axes=axs[0])
    plt.tight_layout()

    # reset branch length to time (in substitution rate units)
    for n in ebola.tree.find_clades():
        n.branch_length=n.clock_length

    # OUTPUT the GTR model
    print(ebola.gtr)

    # plot Skyline
    skyline, confidence = ebola.merger_model.skyline_inferred(gen=50, confidence=2.0)
    skyline_empirical = ebola.merger_model.skyline_empirical(gen=50)
    plt.figure()
    plt.fill_between(skyline.x, confidence[0], confidence[1], color=(0.8, 0.8, 0.8))
    plt.plot(skyline.x, skyline.y, label='maximum likelihood skyline')
    plt.plot(skyline_empirical.x, skyline_empirical.y, label='empirical skyline')
    plt.yscale('log')
    plt.legend()
    plt.ticklabel_format(axis='x',useOffset=False)





