from __future__ import print_function, division
import numpy as np
from Bio import Phylo
from treetime import TreeTime

import matplotlib.pyplot as plt
from matplotlib import cm
try:
    import seaborn as sns
    sns.set_style('whitegrid')
except:
    print ("Seaborn not found. Default style will be used for the plots")


def read_dates(base_name):
    with open(base_name+'.metadata.csv') as date_file:
        dates = {}
        for line in date_file:
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue
    return dates

def format_axes(fig, axs):
    axs[0].set_axis_off()
    axs[1].tick_params(labelsize=14)
    axs[1].set_ylabel("root-to-tip distance", fontsize=16)
    axs[1].set_xlabel("date", fontsize=16)
    fig.tight_layout()


if __name__ == '__main__':

    # load data and parse dates
    plt.ion()
    base_name = 'data/H3N2_NA_allyears_NA.20'
    dates = read_dates(base_name)

    # define a multiplicity for each node (using season here, could be anything)
    seq_multiplicity = {name: 5*(np.cos(2*np.pi*(dates[name]+0.5))+1.2) for name in dates}

    # multiplicity is passed into treetime as a dictionary linking node name to count
    tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', seq_multiplicity=seq_multiplicity,
                          aln = base_name+'.fasta', verbose = 1, dates = dates)

    tt.reroot(root="best")
    fig, axs = plt.subplots(1,2, figsize=(18,9))
    axs[0].set_title("Tree rerooted by treetime", fontsize=18)
    axs[1].set_title("Optimal divergence-time relationship with weighed nodes", fontsize=18)
    Phylo.draw(tt.tree, show_confidence=False, axes=axs[0], label_func=lambda x:x.name.split('|')[0] if x.is_terminal() else "")

    d = np.array([(n.numdate_given, n.dist2root, n.count) for n in tt.tree.get_terminals()])
    mean_count = d[:,2].mean()
    plt.scatter(d[:,0], d[:,1], s=d[:,2]*30/mean_count)
    format_axes(fig, axs)

