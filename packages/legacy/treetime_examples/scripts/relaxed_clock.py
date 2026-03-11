"""
This script illustrates the use of the relaxed clock model in TreeTime.
TreeTime currently only supports a normally distributed clock model.
"""

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

    # load data and parse dates
    plt.ion()
    base_name = '../data/h3n2_na/h3n2_na_20'
    dates = parse_dates(base_name+'.metadata.csv')

    # instantiate treetime
    tt_relaxed = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                        aln = base_name+'.fasta', verbose = 4, dates = dates)

    # this example uses an autocorrelated molecular clock with normal prior and parent-child coupling
    # the parameter slack penalizes rate deviations from the average rate
    # couplings penalize rate changes between parent and child nodes.
    tt_relaxed.run(root='best', relaxed_clock={"slack":5.0, "coupling":1.0}, max_iter=3,
                   resolve_polytomies=True, Tc=0, time_marginal=False)



    ##############
    # PLOTTING
    ##############

    # draw trees inferred with the relaxed model
    fig = plt.figure()
    ax = plt.subplot(111)
    vmin, vmax = 0.5, 1.5 # color branches according to the rate deviation
    for n in tt_relaxed.tree.find_clades():
        if n.up:
            n.color = [int(x*255) for x in cm.cool((min(max(vmin, n.branch_length_interpolator.gamma),vmax)-vmin)/(vmax-vmin))[:3]]
        else:
            n.color = [200,200,200]

    ax.set_title("relaxed clock")
    Phylo.draw(tt_relaxed.tree, axes=ax, show_confidence=False, label_func = lambda x:'')

    # Scatter branch stretch against the rate multiplier of the branch.
    # they are expected to have a positive correlation
    # 1) collect the optimal branch lenght (called mutation_length),
    #    the timetree branch length (clock_length) and the inferred rate multiplier gamma
    branch_lengths = []
    for n in tt_relaxed.tree.find_clades():
        if n.up:
            branch_lengths.append((n.mutation_length, n.clock_length, n.branch_length_interpolator.gamma))


    # 2) plot the difference between optimal and timetree branch length vs the rate multiplier
    plt.figure()
    branch_lengths = np.array(branch_lengths)
    plt.scatter(branch_lengths[:,0]-branch_lengths[:,1], branch_lengths[:,2])
    plt.xlabel("stretch")
    plt.ylabel("rate multiplier")


