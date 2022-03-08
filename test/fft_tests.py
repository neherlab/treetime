import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import pandas as pd

from treetime import TreeTime
from treetime.utils import parse_dates


def get_test_node(tree_list, pattern):
    return [[n for n in tt.tree.find_clades() if pattern in n.name][0] for tt in tree_list]

def get_tree_events(tt):
    tree_events_tt = sorted([(n.time_before_present, n.name, int(n.bad_branch)) for n in tt.tree.find_clades()
                            if not n.bad_branch], key=lambda x:-x[0])
    return tree_events_tt

def compare(new_times_and_names, old_times_and_names):
    new_times_and_names = pd.DataFrame(new_times_and_names)
    old_times_and_names = pd.DataFrame(old_times_and_names)
    old_times_and_names = old_times_and_names.set_index(1)
    new_times_and_names = new_times_and_names.set_index(1)
    old_times_and_names = old_times_and_names.reindex(index=new_times_and_names.index)
    bad_branches = np.float_(old_times_and_names.iloc[:,1]) + 2*np.float_(new_times_and_names.iloc[:,1])
    df = pd.DataFrame(dict(time=np.float_(new_times_and_names.iloc[:,0]),
        difference=(np.float_(new_times_and_names.iloc[:,0])-np.float_(old_times_and_names.iloc[:,0])),
        bad_branch=bad_branches), index=new_times_and_names.index)
    return np.all(new_times_and_names.iloc[:,0]==old_times_and_names.iloc[:,0]), df


if __name__ == '__main__':
    plt.ion()

    ##model parameters for testing
    # choose if should be tested on ebola or h3n2_na dataset and if this script is run
    # on the masterbranch or a branch to be tested
    ebola=True
    master = False

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
    tt_fft = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', use_fft=True,
                    aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)

    for tree in [tt, tt_fft]:
        tree._set_branch_length_mode(seq_kwargs["branch_length_mode"])
        tree.infer_ancestral_sequences(infer_gtr=False, marginal=seq_kwargs["marginal_sequences"])
        tree.prune_short_branches()
        tree.clock_filter(reroot='least-squares', n_iqd=1, plot=False, fixed_clock_rate=tt_kwargs["clock_rate"])
        tree.reroot(root='least-squares', clock_rate=tt_kwargs["clock_rate"])
        tree.infer_ancestral_sequences(**seq_kwargs)
        tree.make_time_tree(clock_rate=tt_kwargs["clock_rate"], time_marginal=tt_kwargs["time_marginal"])
    ##should be no difference at this point unless 'joint' is used for "branch_length_mode"
    tree_events_tt = get_tree_events(tt)
    tree_events_tt_fft = get_tree_events(tt_fft)

    if coal_kwargs['time_marginal']=='assign':
        test_node, test_node_fft = get_test_node([tt, tt_fft], node_pattern)
        while test_node:
            if test_node.name != tt.tree.root.name and not test_node.marginal_pos_LH.is_delta:
                plt.figure()
                plt.plot(test_node.marginal_pos_LH.x, test_node.marginal_pos_LH.prob_relative(test_node.marginal_pos_LH.x),marker='o', linestyle='', ms=2, label='new')
                plt.plot(test_node_fft.marginal_pos_LH.x, test_node_fft.marginal_pos_LH.prob_relative(test_node_fft.marginal_pos_LH.x), marker='o', linestyle='', ms=1, label='fft')
                plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.0001))
                plt.title(test_node.name)
                plt.legend()
                plt.show()
            test_node, test_node_fft= test_node.up, test_node_fft.up

    output_comparison = compare(tree_events_tt_fft, tree_events_tt)
    groups = output_comparison[1].groupby('bad_branch')
    plt.figure()
    for name, group in groups:
        plt.plot(output_comparison[1].time, output_comparison[1].difference, marker='o', linestyle='', ms=1, label=name)
    plt.xlabel("nodes ranging from root at 0 to most recent")
    plt.ylabel("difference time_before_present coalescent branch - master branch")

    large_differences = output_comparison[1][abs(output_comparison[1]['difference']) > abs(np.mean(output_comparison[1].difference)) + np.std(output_comparison[1].difference)]

    for n in list(large_differences.index):
        test_node, test_node_fft = get_test_node([tt, tt_fft], n)
        print([c.name for c in test_node.clades])
        print([c.name for c in test_node_fft.clades])
        print(test_node.up.name + test_node_fft.up.name)
        if test_node.name != tt.tree.root.name and not test_node.marginal_pos_LH.is_delta:
            plt.figure()
            plt.plot(test_node.marginal_pos_LH.x, test_node.marginal_pos_LH.prob_relative(test_node.marginal_pos_LH.x), marker='o', linestyle='', ms=2, label='numerical')
            plt.plot(test_node_fft.marginal_pos_LH.x, test_node_fft.marginal_pos_LH.prob_relative(test_node_fft.marginal_pos_LH.x), marker='o', linestyle='', ms=1, label='fft')
            plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.0001))
            plt.title(test_node.name)
            plt.legend()
            plt.show()
    tt.add_coalescent_model(coal_kwargs ["Tc"])
    tt.make_time_tree(clock_rate=tt_kwargs ["clock_rate"], time_marginal=coal_kwargs ["time_marginal"])
    tree_events_tt_post_coal = get_tree_events(tt)

