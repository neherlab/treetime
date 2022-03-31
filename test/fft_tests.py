## script to compare numerical to fft calculation of convolution integrals,
## will look closer at distributions at nodes where there is a large difference
## between the two different inferred divergence times

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from treetime import TreeTime
from treetime.utils import parse_dates
from treetime.distribution import Distribution


def get_test_nodes(tree_list, pattern):
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

def plot_differences(df, title=None):
    groups = df.groupby('bad_branch')
    fig = plt.figure()
    for name, group in groups:
        plt.plot(df.time, df.difference, marker='o', linestyle='', ms=1, label=name)
    plt.xlabel("nodes ranging from root at 0 to most recent")
    plt.ylabel("difference time_before_present fft - numerical")
    if title:
        plt.title(title)
    plt.show()
    return fig

def compare_dist(node_num, node_fft, dist_node_num, dist_node_fft, dist_name, xlimits=None):
    fig = plt.figure()
    plt.plot(dist_node_num.x, dist_node_num.prob_relative(dist_node_num.x), marker='o', linestyle='', ms=2, label='numerical')
    plt.plot(dist_node_fft.x, dist_node_fft.prob_relative(dist_node_fft.x), marker='o', linestyle='', ms=1, label='fft')
    if xlimits:
        if xlimits!='empty':
            plt.xlim(xlimits)
    else:
        x_min = min(node_num.time_before_present, node_fft.time_before_present)
        x_max = max(node_num.time_before_present, node_fft.time_before_present)
        plt.xlim((x_min-0.0001, x_max+0.0001))
    plt.title(node_num.name + "_" + dist_name)
    plt.legend()
    plt.show()
    return fig


if __name__ == '__main__':
    plt.ion()

    ##model parameters for testing
    # choose if should be tested on ebola or h3n2_na dataset
    ebola=True

    if ebola:
        base_name = '../treetime_examples/data/ebola/ebola'
        clock_rate = 0.0001
    else:
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
        tree.make_time_tree(clock_rate=tt_kwargs["clock_rate"], time_marginal=tt_kwargs["time_marginal"])

    tree_events_tt = get_tree_events(tt)
    tree_events_tt_fft = get_tree_events(tt_fft)

    output_comparison = compare(tree_events_tt_fft, tree_events_tt)
    plot_differences(output_comparison[1])

    large_differences = get_large_differences(output_comparison[1])

    def closer_look_differences():
        ##look closer at the distribution objects where there is a large difference between the two different methods
        for n in list(large_differences.index):
            test_node, test_node_fft = get_test_nodes([tt, tt_fft], n)
            if test_node.name != tt.tree.root.name and not test_node.marginal_pos_LH.is_delta:
                #compare marginal_pos_LH
                compare_dist(test_node, test_node_fft, test_node.marginal_pos_LH, test_node_fft.marginal_pos_LH, 'pos_LH')
                #compare msg_from_parent
                compare_dist(test_node, test_node_fft, test_node.msg_from_parent, test_node_fft.msg_from_parent, 'msg_from_parent')
                #compare subtree_distribution
                compare_dist(test_node, test_node_fft, test_node.subtree_distribution, test_node_fft.subtree_distribution, 'subtree_dist')
                #compare branch_length_interpolator
                compare_dist(test_node, test_node_fft, test_node.branch_length_interpolator, test_node_fft.branch_length_interpolator, 'branch length dist', xlimits='empty')
                #compare marginal_pos_Lx
                xlimits = (test_node.up.time_before_present-0.01,test_node.up.time_before_present+0.01)
                compare_dist(test_node, test_node_fft, test_node.marginal_pos_Lx, test_node_fft.marginal_pos_Lx, 'marginal_pos_Lx', xlimits=xlimits)

                #compare marginal_pos_LH of parent
                print("The parent of numerical is: " + str(test_node.up.name) + " the parent of fft is: " + str(test_node_fft.up.name))
                compare_dist(test_node.up, test_node_fft.up, test_node.up.marginal_pos_LH, test_node_fft.up.marginal_pos_LH, 'pos_LH parent')

                #compare marginal_pos_LH of children
                print("Now looking at children of node of interest")
                for c in test_node.clades:
                    test_node, test_node_fft = get_test_nodes([tt, tt_fft], c.name)
                    print("Time difference is " + str(test_node.time_before_present - test_node_fft.time_before_present ))
                    print("Time of numerical:" + str(test_node.time_before_present) + " time of fft:"+ str(test_node_fft.time_before_present))
                    if test_node.marginal_pos_LH.is_delta:
                        print("delta")
                    else:
                        compare_dist(test_node, test_node_fft, test_node.marginal_pos_LH, test_node_fft.marginal_pos_LH, 'pos_LH')

    closer_look_differences()

    ##now do the same for when the coalescent model is added
    for tree in [tt, tt_fft]:
        tree.add_coalescent_model(coal_kwargs ["Tc"])
        tree.make_time_tree(clock_rate=tt_kwargs ["clock_rate"], time_marginal=coal_kwargs ["time_marginal"])

    tree_events_tt = get_tree_events(tt)
    tree_events_tt_fft = get_tree_events(tt_fft)

    output_comparison = compare(tree_events_tt_fft, tree_events_tt)
    large_differences = get_large_differences(output_comparison[1])
    #plot_differences(output_comparison[1])
    closer_look_differences()

    nodes_to_look_at = ['NODE_0000120', 'NODE_0000098']
    def compare_multiply_functions(nodes_of_interest):
        ## code to look closer at the msgs of a child of a node and behavior of functions
        for n in nodes_of_interest:
            test_node, test_node_fft = get_test_nodes([tt, tt_fft], n)
            print([c.name for c in test_node.clades])
            print([c.name for c in test_node_fft.clades])
            msgs_to_multiply_tt_fft = [child.marginal_pos_Lx for child in test_node_fft.clades
                                            if child.marginal_pos_Lx is not None]
            msgs_to_multiply_tt = [child.marginal_pos_Lx for child in test_node.clades
                                            if child.marginal_pos_Lx is not None]
            subtree_distribution_tt = Distribution.multiply(msgs_to_multiply_tt)
            subtree_distribution_tt_fft = Distribution.multiply(msgs_to_multiply_tt_fft)
            compare_dist(test_node, test_node_fft, subtree_distribution_tt, subtree_distribution_tt_fft, 'multiplied Lx of children')

    # nodes_to_look_at = ['NODE_0000120', 'NODE_0000098']
    # compare_multiply_functions(nodes_to_look_at)