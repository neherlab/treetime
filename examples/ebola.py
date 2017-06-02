'''
this script illustrates the use of treetime to analyze viral sequences.
As an example, we use Ebola virus sequences from the 2014-2015 outbreak in
West Africa. This example contains more than 300 sequences, it will take
a few minutes to run.
'''

from __future__ import print_function, division
from treetime import TreeTime
import numpy as np
from scipy import optimize as sciopt

if __name__ == '__main__':

    # load data and parse dates
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import seaborn as sns
    sns.set_style('white')
    from Bio import Phylo
    plt.ion()
    base_name = 'data/ebola'
    import datetime
    from treetime.utils import numeric_date
    with open(base_name+'.csv') as date_file:
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
    ebola = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                        aln = base_name+'.fasta', verbose = 4, dates = dates)

    # infer an ebola time tree while rerooting and resolving polytomies
    ebola.run(root='best', relaxed_clock=False, max_iter=2,
              resolve_polytomies=True, Tc='skyline', do_marginal="assign")

    # scatter root to tip divergence vs sampling date
    ebola.plot_root_to_tip(add_internal=True)
    t=np.array([2014,2016])
    plt.plot(t, t*ebola.date2dist.clock_rate+ ebola.date2dist.intercept,
             label="y = %1.5f t%1.3f"%(ebola.date2dist.clock_rate, ebola.date2dist.intercept))
    plt.legend(loc=2)

    # rescale branch length to years and plot in axis 0
    from treetime.io import plot_vs_years
    fig, axs = plt.subplots(1,2, sharey=True, figsize=(12,8))
    plot_vs_years(ebola, years=1, ax=axs[1], confidence=(0.05,0.95), label_func = lambda x:"")
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
        if n.up:
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





