import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib import cm

from treetime import TreeTime
from treetime.utils import parse_dates


def get_test_node(tt, pattern):
    return [n for n in tt.tree.find_clades() if pattern in n.name][0]

def get_tree_events(tt):
    tree_events_tt = sorted([(n.time_before_present, n.name) for n in tt.tree.find_clades()
                            if not n.bad_branch], key=lambda x:-x[0])
    return tree_events_tt

def write_to_file(times_and_names):
    with open('master.txt', 'w') as f:
        f.write("time \t name \n")
        for t in times_and_names:
            f.write(str(t[0]))
            f.write("\t")
            f.write(t[1])
            f.write("\n")
    f.close()
    return None

def read_from_file(file_name):
    times_and_names = []
    with open(file_name, 'r') as f:
        for line in f.readlines()[1:]:
            lines = line.split("\n")[0].split("\t")
            times_and_names += [[float(lines[0]), lines[1]]]
    return times_and_names

def compare(file_name, new_times_and_names):
    import pandas as pd
    new_times_and_names = pd.DataFrame(new_times_and_names)
    old_times_and_names = pd.DataFrame(read_from_file(file_name))
    old_times_and_names = old_times_and_names.set_index(1)
    new_times_and_names = new_times_and_names.set_index(1)
    old_times_and_names = old_times_and_names.reindex(index=new_times_and_names.index)
    return np.all(new_times_and_names==old_times_and_names), np.float_(new_times_and_names.iloc[:,0])-np.float_(old_times_and_names.iloc[:,0])

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
                    'time_marginal':False}
    coal_kwargs ={'Tc':10000,
                    'time_marginal':False}

    dates = parse_dates(base_name+'.metadata.csv')
    tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                    aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)

    tt._set_branch_length_mode(seq_kwargs["branch_length_mode"])
    tt.infer_ancestral_sequences(infer_gtr=False, marginal=seq_kwargs["marginal_sequences"])
    tt.prune_short_branches()
    tt.clock_filter(reroot='least-squares', n_iqd=1, plot=False, fixed_clock_rate=tt_kwargs["clock_rate"])
    tt.reroot(root='least-squares', clock_rate=tt_kwargs["clock_rate"])
    tt.infer_ancestral_sequences(**seq_kwargs)
    tt.make_time_tree(clock_rate=tt_kwargs["clock_rate"], time_marginal=tt_kwargs["time_marginal"])
    ##should be no difference at this point unless 'joint' is used for "branch_length_mode"
    # if master:
    #     write_to_file(tree_events_tt)
    # else:
    #     output_comparison = compare("../../TreeTimeMaster/treetime/master.txt", tree_events_tt)

    tt.add_coalescent_model(coal_kwargs ["Tc"])
    tt.make_time_tree(clock_rate=tt_kwargs ["clock_rate"], time_marginal=coal_kwargs ["time_marginal"])
    tree_events_tt_post_coal = get_tree_events(tt)

    if master:
        write_to_file(tree_events_tt_post_coal)
    else:
        output_comparison = compare("../../TreeTimeMaster/treetime/master.txt", tree_events_tt_post_coal)
        plt.figure()
        plt.plot(output_comparison[1], 'o')
        plt.xlabel("nodes ranging from root at 0 to most recent")
        plt.ylabel("difference time_before_present coalescent branch - master branch")
        if ebola:
            plt.savefig("time_before_present-differences-ebola.png")
        else:
            plt.savefig("time_before_present-differences-h3n2_na.png")
    Phylo.draw(tt.tree, label_func=lambda x:"")

    if coal_kwargs['time_marginal']=='assign':
        test_node = get_test_node(tt, node_pattern)
        while test_node:
            if test_node.name != tt.tree.root.name:
                plt.figure()
                t = np.linspace(test_node.time_before_present-0.0001,test_node.time_before_present+0.0001,1000)
                plt.plot(t, test_node.marginal_pos_LH.prob_relative(t), label='new', ls='--')
                plt.title(test_node.name)
                plt.show()
            test_node= test_node.up