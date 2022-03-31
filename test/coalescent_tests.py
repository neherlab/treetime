## This script is written for running two different versions of treetime (in different environments)
## and comparing their output, in this case: inferred divergence times of nodes.
## Remember to set if this script is currently being run on the master branch or
## the new branch/version that should be tested

import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import pandas as pd

from treetime import TreeTime
from treetime.utils import parse_dates
from fft_tests import get_tree_events, get_test_nodes, get_large_differences, compare


def write_to_file(times_and_names, name = 'master.txt'):
    with open(name, 'w') as f:
        f.write("time \t name \t bad_branch\n")
        for t in times_and_names:
            f.write(str(t[0]))
            f.write("\t")
            f.write(t[1])
            f.write("\t")
            f.write(str(t[2]))
            f.write("\n")
    f.close()
    return None

def read_from_file(file_name):
    times_and_names = []
    with open(file_name, 'r') as f:
        for line in f.readlines()[1:]:
            lines = line.split("\n")[0].split("\t")
            times_and_names += [[float(lines[0]), lines[1], lines[2]]]
    return times_and_names



if __name__ == '__main__':
    plt.ion()

    ebola=True
    master = False ##should be True for master branch
    location_master = '../../TreeTimeMaster/treetime/' ##only needed for test branch

    if ebola:
        node_pattern = 'EM_004555'
        base_name = '../treetime_examples/data/ebola/ebola'
        clock_rate = 0.0001
    else:
        node_pattern = 'Indiana'
        base_name = '../treetime_examples/data/h3n2_na/h3n2_na_20'
        clock_rate = 0.0028

    seq_kwargs = {"marginal_sequences":True,
                    "branch_length_mode": 'input',
                    "sample_from_profile":"root",
                    "reconstruct_tip_states":False}
    tt_kwargs = {'clock_rate':clock_rate,
                    'time_marginal':'assign'}
    coal_kwargs ={'Tc':10000,
                    'time_marginal':'assign'}

    dates = parse_dates(base_name+'.metadata.csv')
    tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', use_fft=False,
                    aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)

    tt._set_branch_length_mode(seq_kwargs["branch_length_mode"])
    tt.infer_ancestral_sequences(infer_gtr=False, marginal=seq_kwargs["marginal_sequences"])
    tt.prune_short_branches()
    tt.clock_filter(reroot='least-squares', n_iqd=1, plot=False, fixed_clock_rate=tt_kwargs["clock_rate"])
    tt.reroot(root='least-squares', clock_rate=tt_kwargs["clock_rate"])
    tt.infer_ancestral_sequences(**seq_kwargs)
    tt.make_time_tree(clock_rate=tt_kwargs["clock_rate"], time_marginal=tt_kwargs["time_marginal"])
    ##should be no difference at this point unless 'joint' is used for "branch_length_mode"

    tree_events_tt = get_tree_events(tt)
    if master:
        write_to_file(tree_events_tt)
    else:
        output_comparison = compare(tree_events_tt, read_from_file(location_master + "master.txt"))
        large_differences = get_large_differences(output_comparison[1])
        ##plot differences for non bad-branches
        plt.figure()
        plt.plot(output_comparison[1][output_comparison[1]['bad_branch']==0].time, output_comparison[1][output_comparison[1]['bad_branch']==0].difference, 'o')

    # add coalescent model, assume there will be some changes after this point
    tt.add_coalescent_model(coal_kwargs ["Tc"])
    tt.make_time_tree(clock_rate=tt_kwargs ["clock_rate"], time_marginal=coal_kwargs ["time_marginal"])
    tree_events_tt_post_coal = get_tree_events(tt)

    if master:
        write_to_file(tree_events_tt_post_coal)
    else:
        write_to_file(tree_events_tt, "fft_branch"+  ".txt")
        output_comparison = compare(tree_events_tt_post_coal, read_from_file(location_master + "master.txt"))
        ##plot differences and color according to bad-branch label, save figure
        plt.figure()
        groups = output_comparison[1].groupby('bad_branch')
        color = ['red', 'black', 'yellow', 'green']
        for name, group in groups:
            plt.plot(output_comparison[1].time, output_comparison[1].difference, marker='o', color=color[int(name)], linestyle='', ms=1, label=name)
        plt.xlabel("nodes ranging from root at 0 to most recent")
        plt.ylabel("difference time_before_present coalescent branch - master branch")
        plt.legend()
        if ebola:
            plt.savefig("time_before_present-differences-ebola.png")
        else:
            plt.savefig("time_before_present-differences-h3n2_na.png")
        large_differences = get_large_differences(output_comparison[1])
        ##plot differences for non bad-branches
        plt.figure()
        plt.plot(output_comparison[1][output_comparison[1]['bad_branch']==0].time, output_comparison[1][output_comparison[1]['bad_branch']==0].difference, 'o')

    ##draw tree for optical comparison
    Phylo.draw(tt.tree, label_func=lambda x:"")

    ##plot LH distributions for optical comparison
    if coal_kwargs['time_marginal']=='assign':
        test_node = get_test_nodes([tt], node_pattern)[0]
        while test_node:
            if test_node.name != tt.tree.root.name:
                plt.figure()
                t = np.linspace(test_node.time_before_present-0.0001,test_node.time_before_present+0.0001,1000)
                plt.plot(t, test_node.marginal_pos_LH.prob_relative(t), 'o-', label='new')
                plt.title(test_node.name)
                plt.show()
            test_node= test_node.up