import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
import pandas as pd

from treetime import TreeTime
from treetime.utils import parse_dates


def get_test_node(tt, pattern):
    return [n for n in tt.tree.find_clades() if pattern in n.name][0]

def get_tree_events(tt):
    tree_events_tt = sorted([(n.time_before_present, n.name, int(n.bad_branch)) for n in tt.tree.find_clades()],
                         key=lambda x:-x[0])
    return tree_events_tt

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

def compare(file_name, new_times_and_names):
    new_times_and_names = pd.DataFrame(new_times_and_names)
    old_times_and_names = pd.DataFrame(read_from_file(file_name))
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
        output_comparison = compare("../../TreeTimeMaster/treetime/master.txt", tree_events_tt)
        large_differences = output_comparison[1][abs(output_comparison[1]['difference']) > abs(np.mean(output_comparison[1].difference)) + 2*np.std(output_comparison[1].difference)]
        plt.figure()
        plt.plot(output_comparison[1][output_comparison[1]['bad_branch']==0].time, output_comparison[1][output_comparison[1]['bad_branch']==0].difference, 'o')


    tt.add_coalescent_model(coal_kwargs ["Tc"])
    tt.make_time_tree(clock_rate=tt_kwargs ["clock_rate"], time_marginal=coal_kwargs ["time_marginal"])
    tree_events_tt_post_coal = get_tree_events(tt)

    if master:
        write_to_file(tree_events_tt_post_coal)
    else:
        write_to_file(tree_events_tt, "fft_branch"+  ".txt")
        output_comparison = compare("../../TreeTimeMaster/treetime/master.txt", tree_events_tt_post_coal)
        groups = output_comparison[1].groupby('bad_branch')
        plt.figure()
        color = ['red', 'black', 'yellow', 'green']
        for name, group in groups:
            print(color[int(name)])
            plt.plot(output_comparison[1].time, output_comparison[1].difference, marker='o', color=color[int(name)], linestyle='', ms=1, label=name)
        plt.xlabel("nodes ranging from root at 0 to most recent")
        plt.ylabel("difference time_before_present coalescent branch - master branch")
        plt.legend()
        if ebola:
            plt.savefig("time_before_present-differences-ebola.png")
        else:
            plt.savefig("time_before_present-differences-h3n2_na.png")
        large_differences = output_comparison[1][abs(output_comparison[1]['difference']) > abs(np.mean(output_comparison[1].difference)) + 2*np.std(output_comparison[1].difference)]
        plt.figure()
        plt.plot(output_comparison[1][output_comparison[1]['bad_branch']==0].time, output_comparison[1][output_comparison[1]['bad_branch']==0].difference, 'o')

    Phylo.draw(tt.tree, label_func=lambda x:"")

    if coal_kwargs['time_marginal']=='assign':
        test_node = get_test_node(tt, node_pattern)
        while test_node:
            if test_node.name != tt.tree.root.name:
                plt.figure()
                t = np.linspace(test_node.time_before_present-0.0001,test_node.time_before_present+0.0001,1000)
                plt.plot(t, test_node.marginal_pos_LH.prob_relative(t), 'o-', label='new')
                plt.title(test_node.name)
                plt.show()
            test_node= test_node.up