from __future__ import print_function, division
from treetime import TreeTime
import numpy as np
from scipy import optimize as sciopt

if __name__ == '__main__':

    # load data
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import seaborn as sns
    sns.set_style('whitegrid')
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

    tt_relaxed = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                        aln = base_name+'.fasta', verbose = 4, dates = dates, debug=True)

    # this example uses an autocorrelated molecular clock with normal prior and parent-child coupling
    tt_relaxed.run(root='best', relaxed_clock={"slack":5.0, "coupling":1.0}, max_iter=1,
               resolve_polytomies=True, Tc=0, do_marginal=False)

    tt_strict = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                        aln = base_name+'.fasta', verbose = 4, dates = dates, debug=True)

    # this example uses an autocorrelated molecular clock with location prior
    tt_strict.run(root='best', relaxed_clock=False, max_iter=1,
               resolve_polytomies=True, Tc=0, do_marginal=False)

    # draw trees inferred with relaxed and strict model
    fig, axs = plt.subplots(1,2, sharex=True, sharey=True, figsize=(12,8))
    vmin=0
    vmax=2.0
    for n in tt_relaxed.tree.find_clades():
        if n.up:
            n.color = [int(x*255) for x in cm.cool((min(max(vmin, n.branch_length_interpolator.gamma),vmax)-vmin)/(vmax-vmin))[:3]]
        else:
            n.color = [200,200,200]
    for n in tt_strict.tree.find_clades():
        n.color = [200,200,200]

    axs[0].set_title("relaxed clock")
    Phylo.draw(tt_strict.tree, axes=axs[1], show_confidence=False, label_func = lambda x:'')
    axs[1].set_title("strict clock")
    Phylo.draw(tt_relaxed.tree, axes=axs[0], show_confidence=False, label_func = lambda x:"")

    # Scatter branch stretch against the rate multiplier of the branch.
    # this is expected to have a positive relationship
    branch_lengths = []
    for n in tt_relaxed.tree.find_clades():
        if n.up:
            branch_lengths.append((n.mutation_length, n.clock_length, n.branch_length_interpolator.gamma))

    branch_lengths = np.array(branch_lengths)

    plt.figure()
    plt.scatter(branch_lengths[:,0]-branch_lengths[:,1], branch_lengths[:,2])
    plt.xlabel("stretch")
    plt.ylabel("rate multiplier")


