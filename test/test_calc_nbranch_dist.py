import matplotlib.pyplot as plt
import numpy as np
from math import floor, ceil

from treetime import TreeTime as TreeTime
from treetime.utils import parse_dates
import treetime.config as ttconf
from scipy.interpolate import interp1d


def undo_merger(test_node, tt, t):
    return test_node.branch_length_interpolator.prob(t)*np.exp(tt.merger_model.cost(test_node.time_before_present, t))

def get_test_nodes(tree_list, pattern):
    return [[n for n in tt.tree.find_clades() if pattern in n.name][0] for tt in tree_list]


if __name__ == '__main__':
    plt.ion()
    base_name = '../treetime_examples/data/h3n2_na/h3n2_na_20'
    #base_name = '../treetime_examples/data/ebola/ebola'

    dates = parse_dates(base_name+'.metadata.csv')

    tt_old= TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                            aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)
    tt_smooth= TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                            aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)

    fixed_clock_rate = 0.0028
    for tt in [tt_old, tt_smooth]:
        tt.reroot(root='least-squares', clock_rate=fixed_clock_rate)
        tt.infer_ancestral_sequences(infer_gtr=False, marginal=False)
        tt.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=False)
        #tt._ml_t_marginal(assign_dates=True)
        tt.add_coalescent_model(0.01)

    tt_smooth.merger_model.calc_branch_count_dist()
    tt_smooth.merger_model.set_Tc(0.01)

    node_pattern = 'Indiana'
    test_nodes = get_test_nodes([tt_old, tt_smooth], node_pattern)
    while test_nodes[0] is not None:
        if test_nodes[0].name != tt.tree.root.name:
            plt.figure()
            t = np.linspace(0,4*fixed_clock_rate,1000)
            plt.plot(t, undo_merger(test_nodes[0], tt_old, t), label='old', ls='-')
            plt.plot(t, undo_merger(test_nodes[1], tt_smooth, t), label='old', ls='-')
            plt.savefig("new.png")
        test_nodes =[test_node.up for test_node in test_nodes]

    plt.figure()
    plt.plot(tt_old.merger_model.nbranches.x, tt_old.merger_model.nbranches.y, label="old nbranches function")
    plt.plot(tt_smooth.merger_model.nbranches.x, tt_smooth.merger_model.nbranches.y, label="new nbranches function")
    plt.legend()
    plt.xlim((0,0.08))