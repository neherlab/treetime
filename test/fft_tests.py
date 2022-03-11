import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import pandas as pd

from treetime import TreeTime
from treetime.utils import parse_dates
from treetime.distribution import Distribution


def get_test_node(tree_list, pattern):
    return [[n for n in tt.tree.find_clades() if pattern in n.name][0] for tt in tree_list]

def get_tree_events(tt):
    tree_events_tt = sorted([(n.time_before_present, n.name, int(n.bad_branch)) for n in tt.tree.find_clades()], key=lambda x:-x[0])
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

def get_large_differences(df):
    large_differences = df[abs(df.difference) > abs(np.mean(df.difference)) + 2*np.std(df.difference)]
    return large_differences




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
    tt_fft = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', use_fft=True, precision_fft=200,
                    aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)

    for tree in [tt, tt_fft]:
        tree._set_branch_length_mode(seq_kwargs["branch_length_mode"])
        tree.infer_ancestral_sequences(infer_gtr=False, marginal=seq_kwargs["marginal_sequences"])
        tree.prune_short_branches()
        tree.clock_filter(reroot='least-squares', n_iqd=1, plot=False, fixed_clock_rate=tt_kwargs["clock_rate"])
        tree.reroot(root='least-squares', clock_rate=tt_kwargs["clock_rate"])
        tree.infer_ancestral_sequences(**seq_kwargs)

    tt.make_time_tree(clock_rate=tt_kwargs["clock_rate"], time_marginal=tt_kwargs["time_marginal"], divide=False)
    tt_fft.make_time_tree(clock_rate=tt_kwargs["clock_rate"], time_marginal=tt_kwargs["time_marginal"], divide=False)
    ##should be no difference at this point unless 'joint' is used for "branch_length_mode"
    tree_events_tt = get_tree_events(tt)
    tree_events_tt_fft = get_tree_events(tt_fft)

    # if coal_kwargs['time_marginal']=='assign':
    #     test_node, test_node_fft = get_test_node([tt, tt_fft], node_pattern)
    #     while test_node:
    #         if test_node.name != tt.tree.root.name and not test_node.marginal_pos_LH.is_delta:
    #             plt.figure()
    #             plt.plot(test_node.marginal_pos_LH.x, test_node.marginal_pos_LH.prob_relative(test_node.marginal_pos_LH.x),marker='o', linestyle='', ms=2, label='new')
    #             plt.plot(test_node_fft.marginal_pos_LH.x, test_node_fft.marginal_pos_LH.prob_relative(test_node_fft.marginal_pos_LH.x), marker='o', linestyle='', ms=1, label='fft')
    #             plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.0001))
    #             plt.title(test_node.name)
    #             plt.legend()
    #             plt.show()
    #         test_node, test_node_fft= test_node.up, test_node_fft.up

    output_comparison = compare(tree_events_tt_fft, tree_events_tt)
    groups = output_comparison[1].groupby('bad_branch')
    plt.figure()
    for name, group in groups:
        plt.plot(output_comparison[1].time, output_comparison[1].difference, marker='o', linestyle='', ms=1, label=name)
    plt.xlabel("nodes ranging from root at 0 to most recent")
    plt.ylabel("difference time_before_present coalescent branch - master branch")

    large_differences = get_large_differences(output_comparison[1])

    for n in list(large_differences.index):
        test_node, test_node_fft = get_test_node([tt, tt_fft], n)
        print([c.name for c in test_node.clades])
        print([c.name for c in test_node_fft.clades])
        if test_node.name != tt.tree.root.name and not test_node.marginal_pos_LH.is_delta:
            plt.figure()
            plt.plot(test_node.marginal_pos_LH.x, test_node.marginal_pos_LH.prob_relative(test_node.marginal_pos_LH.x), marker='o', linestyle='', ms=2, label='numerical')
            plt.plot(test_node_fft.marginal_pos_LH.x, test_node_fft.marginal_pos_LH.prob_relative(test_node_fft.marginal_pos_LH.x), marker='o', linestyle='', ms=1, label='fft')
            plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.0001))
            plt.title(test_node.name + "LH")
            plt.legend()
            plt.show()
            plt.figure()
            plt.plot(test_node.msg_from_parent.x, test_node.msg_from_parent.prob_relative(test_node.msg_from_parent.x), marker='o', linestyle='', ms=2, label='numerical')
            plt.plot(test_node_fft.msg_from_parent.x, test_node_fft.msg_from_parent.prob_relative(test_node_fft.msg_from_parent.x), marker='o', linestyle='', ms=1, label='fft')
            plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.0001))
            plt.title(test_node.name + "msg parent")
            plt.legend()
            plt.show()
            plt.figure()
            plt.plot(test_node.subtree_distribution.x, test_node.subtree_distribution.prob_relative(test_node.subtree_distribution.x), marker='o', linestyle='', ms=2, label='numerical')
            plt.plot(test_node_fft.subtree_distribution.x, test_node_fft.subtree_distribution.prob_relative(test_node_fft.subtree_distribution.x), marker='o', linestyle='', ms=1, label='fft')
            plt.xlim((test_node.time_before_present-0.01,test_node.time_before_present+0.01))
            plt.title(test_node.name + "subtree_dist")
            plt.legend()
            plt.show()
            plt.figure()
            plt.plot(test_node.branch_length_interpolator.x, test_node.branch_length_interpolator.prob_relative(test_node.branch_length_interpolator.x), marker='o', linestyle='', ms=2, label='numerical')
            plt.plot(test_node_fft.branch_length_interpolator.x, test_node_fft.branch_length_interpolator.prob_relative(test_node_fft.branch_length_interpolator.x), marker='o', linestyle='', ms=1, label='fft')
            plt.title(test_node.name + "bracnh length")
            plt.legend()
            plt.show()
            plt.figure()
            plt.plot(test_node.marginal_pos_Lx.x, test_node.marginal_pos_Lx.prob_relative(test_node.marginal_pos_Lx.x), marker='o', linestyle='', ms=2, label='numerical')
            plt.plot(test_node_fft.marginal_pos_Lx.x, test_node_fft.marginal_pos_Lx.prob_relative(test_node_fft.marginal_pos_Lx.x), marker='o', linestyle='', ms=1, label='fft')
            plt.xlim((test_node.up.time_before_present-0.01,test_node.up.time_before_present+0.01))
            plt.title(test_node.name + "marginal Lx")
            plt.legend()
            plt.show()
            print(test_node.up.name + test_node_fft.up.name)
            plt.figure()
            plt.plot(test_node.up.marginal_pos_LH.x, test_node.up.marginal_pos_LH.prob_relative(test_node.up.marginal_pos_LH.x), marker='o', linestyle='', ms=2, label='numerical')
            plt.plot(test_node_fft.up.marginal_pos_LH.x, test_node_fft.up.marginal_pos_LH.prob_relative(test_node_fft.up.marginal_pos_LH.x), marker='o', linestyle='', ms=1, label='fft')
            plt.xlim((test_node.up.time_before_present-0.0001,test_node.up.time_before_present+0.0001))
            plt.title(test_node.up.name)
            plt.legend()
            plt.show()
            msg_parent_to_node = Distribution.divide(test_node.up.marginal_pos_LH, test_node.marginal_pos_Lx)
            msg_parent_to_node_fft = Distribution.divide(test_node_fft.up.marginal_pos_LH, test_node_fft.marginal_pos_Lx)
            plt.figure()
            plt.plot(msg_parent_to_node.x, msg_parent_to_node.prob_relative(msg_parent_to_node.x), marker='o', linestyle='', ms=2, label='numerical')
            plt.plot(msg_parent_to_node_fft.x, msg_parent_to_node_fft.prob_relative(msg_parent_to_node_fft.x), marker='o', linestyle='', ms=1, label='fft')
            plt.title(test_node.up.name)
            plt.legend()
            plt.show()
            for c in test_node.clades:
                test_node, test_node_fft = get_test_node([tt, tt_fft], c.name)
                print("Time difference is " + str(test_node.time_before_present - test_node_fft.time_before_present ))
                print("Time of numerical:" + str(test_node.time_before_present) + " time of fft:"+ str(test_node_fft.time_before_present))
                if test_node.marginal_pos_LH.is_delta:
                    print("delta")
                else:
                    plt.figure()
                    plt.plot(test_node.marginal_pos_LH.x, test_node.marginal_pos_LH.prob_relative(test_node.marginal_pos_LH.x), marker='o', linestyle='', ms=2, label='numerical')
                    plt.plot(test_node_fft.marginal_pos_LH.x, test_node_fft.marginal_pos_LH.prob_relative(test_node_fft.marginal_pos_LH.x), marker='o', linestyle='', ms=1, label='fft')
                    plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.001))
                    plt.title(test_node.name)
                    plt.legend()
                    plt.show()

    tt.add_coalescent_model(coal_kwargs ["Tc"])
    tt.make_time_tree(clock_rate=tt_kwargs ["clock_rate"], time_marginal=coal_kwargs ["time_marginal"])
    tt_fft.add_coalescent_model(coal_kwargs ["Tc"])
    tt_fft.make_time_tree(clock_rate=tt_kwargs ["clock_rate"], time_marginal=coal_kwargs ["time_marginal"], divide=False)
    tree_events_tt_post_coal = get_tree_events(tt)

    ## code to look closer at the msgs of a child of a node and behavior of functions
    #nodes_to_look_at = ['NODE_0000120', 'NODE_0000098']
    # for n in nodes_to_look_at:
    #     test_node, test_node_fft = get_test_node([tt, tt_fft], n)
    #     print([c.name for c in test_node.clades])
    #     print([c.name for c in test_node_fft.clades])
    #     msgs_to_multiply_tt_fft = [child.marginal_pos_Lx for child in test_node_fft.clades
    #                                     if child.marginal_pos_Lx is not None]
    #     msgs_to_multiply_tt = [child.marginal_pos_Lx for child in test_node.clades
    #                                     if child.marginal_pos_Lx is not None]
    #     subtree_distribution_tt = Distribution.multiply(msgs_to_multiply_tt)
    #     subtree_distribution_tt_fft = Distribution.multiply(msgs_to_multiply_tt_fft)
    #     plt.figure()
    #     plt.plot(subtree_distribution_tt_fft.x, subtree_distribution_tt_fft.prob_relative(subtree_distribution_tt_fft.x), marker='o', linestyle='', ms=2, label='fft')
    #     plt.plot(subtree_distribution_tt.x, subtree_distribution_tt.prob_relative(subtree_distribution_tt.x), marker='o', linestyle='', ms=1, label='numerical')
    #     plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.001))
    #     plt.title(test_node.name + "child LH")
    #     plt.legend()
    #     plt.show()
    #     if test_node.name != tt.tree.root.name and not test_node.marginal_pos_LH.is_delta:
    #         for c in test_node.clades:
    #             test_node, test_node_fft = get_test_node([tt, tt_fft], c.name)
    #             plt.figure()
    #             plt.plot(test_node.marginal_pos_Lx.x, test_node.marginal_pos_Lx.prob_relative(test_node.marginal_pos_Lx.x), marker='o', linestyle='', ms=2, label='numerical')
    #             plt.plot(test_node_fft.marginal_pos_Lx.x, test_node_fft.marginal_pos_Lx.prob_relative(test_node_fft.marginal_pos_Lx.x), marker='o', linestyle='', ms=1, label='fft')
    #             plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.001))
    #             plt.title(test_node.name + "child Lx")
    #             plt.legend()
    #             plt.show()