import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib import cm

from treetime import TreeTime
from treetime_fft import TreeTime as TreeTimeFFT
from treetime.utils import parse_dates


def get_test_nodes(tt1, tt2, pattern):
    return [n for n in tt1.tree.find_clades() if pattern in n.name][0], [n for n in tt2.tree.find_clades() if pattern in n.name][0]


if __name__ == '__main__':

    # load data and parse dates
    plt.ion()
    base_name = '../treetime_examples/data/h3n2_na/h3n2_na_20'

    dates = parse_dates(base_name+'.metadata.csv')
    tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', use_fft=False,
                  aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3)
    tt_fft = TreeTimeFFT(gtr='Jukes-Cantor', tree = base_name+'.nwk', use_fft=False,
                  aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3)

    tt.min_width = 0.00001
    tt_fft.min_width = 0.00001
    tt.debug = True
    tt_fft.debug = True

    fixed_clock_rate = 0.0028
    tt.reroot(root='least-squares', clock_rate=fixed_clock_rate)
    tt_fft.reroot(root='least-squares', clock_rate=fixed_clock_rate)
    tt.infer_ancestral_sequences(infer_gtr=False, marginal=False)
    tt_fft.infer_ancestral_sequences(infer_gtr=False, marginal=False)

    tt.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=False)
    tt_fft.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=False)

    node_pattern = 'NODE_0000003'
    test_node, test_node_fft = get_test_nodes(tt, tt_fft, node_pattern)

    plt.figure()
    t = np.linspace(0,4*fixed_clock_rate,1000)
    plt.plot(t, test_node.branch_length_interpolator.prob_relative(t), label='old', ls='-')
    plt.plot(t, test_node_fft.branch_length_interpolator.prob_relative(t), label='new', ls='--')

    Tc=0.01
    tt.add_coalescent_model(Tc)
    tt_fft.add_coalescent_model(Tc)

    plt.figure()
    plt.plot(t, tt.merger_model.cost(test_node.time_before_present, t))
    plt.plot(t, tt_fft.merger_model.cost(test_node.time_before_present, t))

    plt.figure()
    plt.plot(t, tt.merger_model.integral_merger_rate(t+0.02))
    plt.plot(t, tt_fft.merger_model.integral_merger_rate(t+0.02))

    plt.figure()
    t = np.linspace(0,4*fixed_clock_rate,1000)
    plt.plot(t, test_node.branch_length_interpolator.prob_relative(t), label='old', ls='-')
    plt.plot(t, test_node_fft.branch_length_interpolator.prob_relative(t), label='new', ls='--')

    plt.figure()
    t = np.linspace(0,4*fixed_clock_rate,1000)
    def undo_merger(test_node, t):
        return test_node.branch_length_interpolator.prob(t)*np.exp(tt.merger_model.cost(test_node.time_before_present, t))

    plt.plot(t, undo_merger(test_node, t), label='old', ls='-')
    plt.plot(t, test_node_fft.branch_length_interpolator.prob(t), label='new', ls='--')


    tt.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=True)
    tt_fft.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=True)


    def undo_merger_H(test_node, t):
        res = test_node.subtree_distribution.prob_relative(t)*np.exp(-tt_fft.merger_model.integral_merger_rate(t))
        return res/res.max()
    while test_node:
        t_node = np.linspace(test_node.time_before_present-0.005,test_node.time_before_present+0.005,1000)
        plt.figure(f"subtree distribution {test_node.name}, mult: {len(test_node.clades)}")

        plt.plot(t_node, test_node.subtree_distribution.prob_relative(t_node), label='old', ls='-')
        plt.plot(t_node, test_node_fft.subtree_distribution.prob_relative(t_node), label='new', ls='-')
        plt.plot(t_node, undo_merger_H(test_node_fft, t_node), label='new_corrected', ls='--')
        plt.legend()

        test_node, test_node_fft = test_node.up, test_node_fft.up

    test_node = tt.tree.root
    test_node_fft = tt_fft.tree.root
    t_node = np.linspace(test_node.time_before_present-0.005,test_node.time_before_present+0.005,1000)

    plt.figure("marginal_pos")
    plt.plot(t_node, test_node.marginal_pos_LH.prob_relative(t_node), label='old', ls='-')
    plt.plot(t_node, test_node_fft.marginal_pos_LH.prob_relative(t_node), label='new', ls='--')

    plt.figure("marginal_pos_corrected")
    def undo_merger_P(test_node, t):
        res = test_node.marginal_pos_LH.prob_relative(t)*np.exp(-tt_fft.merger_model.integral_merger_rate(t))
        return res/res.max()

    plt.plot(t_node, test_node.marginal_pos_LH.prob_relative(t_node), label='old', ls='-')
    plt.plot(t_node, test_node_fft.marginal_pos_LH.prob_relative(t_node), label='new', ls='--')
    plt.plot(t_node, undo_merger_P(test_node_fft, t_node), label='new', ls=':')

