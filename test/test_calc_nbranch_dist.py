import matplotlib.pyplot as plt
import numpy as np
from math import floor, ceil

from treetime import TreeTime as TreeTime
from treetime.utils import parse_dates
import treetime.config as ttconf
from scipy.interpolate import interp1d
from Bio import Phylo


def undo_merger(test_node, tt, t):
    return test_node.branch_length_interpolator.prob(t)*np.exp(tt.merger_model.cost(test_node.time_before_present, t))

def get_test_nodes(tree_list, pattern):
    return [[n for n in tt.tree.find_clades() if pattern in n.name][0] for tt in tree_list]


if __name__ == '__main__':
    plt.ion()

    ebola=True
    if ebola:
        base_name = '../treetime_examples/data/ebola/ebola'
    else:
        base_name = '../treetime_examples/data/h3n2_na/h3n2_na_20'

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

    ## set the branch_count function using the smooth approach and time differences
    import time
    start = time.process_time()
    tt_smooth.add_coalescent_model(Tc=0.001, discrete_nbranches=False)
    print("Time for smooth nbranches:" + str(time.process_time()-start))
    start = time.process_time()
    tt_old.add_coalescent_model(Tc=0.001, discrete_nbranches=True)
    print("Time for discrete nbranches:" + str(time.process_time()-start))
    if ebola:
        node_pattern= 'V517'
    else:
        node_pattern = 'Indiana'

    ##if still checking against old function
    tt_old.merger_model.calc_branch_count()

    ## Plot differences in the nbranches interp1d object
    plt.figure()
    plt.plot(tt_old.merger_model.nbranches.x/tt_old.date2dist.clock_rate, tt_old.merger_model.nbranches.y, label="old nbranches function")
    plt.plot(tt_smooth.merger_model.nbranches.x/tt_smooth.date2dist.clock_rate, tt_smooth.merger_model.nbranches.y, label="new smooth nbranches function")
    plt.xlabel("time before present")
    plt.ylabel("nbranches")
    plt.legend()
    plt.xlim((0,40))

    x_old = tt_smooth.merger_model.nbranches.x/tt_smooth.date2dist.clock_rate
    y_old = tt_smooth.merger_model.nbranches.y
    tt_smooth.merger_model.calc_branch_count_dist(discrete=True)
    ## Plot differences in the nbranches interp1d object
    plt.figure()
    plt.plot(tt_old.merger_model.nbranches.x/tt_old.date2dist.clock_rate, tt_old.merger_model.nbranches.y, label="old nbranches function")
    plt.plot(x_old, y_old, label="new smooth nbranches function")
    plt.plot(tt_smooth.merger_model.nbranches.x/tt_smooth.date2dist.clock_rate, tt_smooth.merger_model.nbranches.y, label="new discrete nbranches function")
    plt.xlabel("time before present")
    plt.ylabel("nbranches")
    plt.legend()
    plt.xlim((0,40))

    ##check that this is really the same as the previously calculated calc_branch_count:
    print(np.all(tt_smooth.merger_model.observed_nbranches.y==tt_old.merger_model.nbranches.y))

    ## Plot effects on branch length distribution and cost function of coalescent
    test_nodes = get_test_nodes([tt_old, tt_smooth], node_pattern)
    while test_nodes[0] is not None:
        if test_nodes[0].name != tt.tree.root.name:
            plt.figure()
            t = np.linspace(0,4*fixed_clock_rate,1000)
            plt.plot(t, undo_merger(test_nodes[0], tt_old, t), label='old', ls='-')
            plt.plot(t, undo_merger(test_nodes[1], tt_smooth, t), label='smooth', ls='-')
            plt.legend()
        test_nodes =[test_node.up for test_node in test_nodes]


    ## Plot effects on final constructed time trees
    for tt in [tt_old, tt_smooth]:
        tt.make_time_tree(clock_rate=fixed_clock_rate, time_marginal=True)

    fig, axs = plt.subplots(1,2, sharey=True, figsize=(12,8))
    Phylo.draw(tt_old.tree, label_func=lambda x:"", axes=axs[0])
    axs[0].set_title("old time tree")
    Phylo.draw(tt_smooth.tree, label_func=lambda x:"", axes=axs[1])
    axs[1].set_title("smooth time tree")

    ## Plot effects on final LH function
    test_nodes = get_test_nodes([tt_old, tt_smooth], node_pattern)
    next=0
    fig, axs = plt.subplots(2,2, sharey=True, figsize=(12,8))
    fig.suptitle("Effect of smooth nbranches on LH distribution")
    while test_nodes[0] is not None:
        if test_nodes[0].name != tt.tree.root.name:
            if next==0:
                i, j=0, 0
            elif next==1:
                i, j=1, 0
            elif next==2:
                i, j=0, 1
            elif next==3:
                i, j=1, 1
            t = np.linspace(test_nodes[0].time_before_present-0.005,test_nodes[0].time_before_present+0.003,2000)
            axs[i,j].plot(t, test_nodes[0].marginal_pos_LH.prob_relative(t),  label='old', ls='-')
            axs[i,j].plot(t, test_nodes[1].marginal_pos_LH.prob_relative(t), label='smooth', ls='-')
            axs[i,j].set_xlabel("time before present")
            axs[i,j].set_ylabel("marginal LH")
            axs[i,j].legend()
            next += 1
            if next==4:
                next=0
                fig, axs = plt.subplots(2,2, sharey=True, figsize=(12,8))
                fig.suptitle("Effect of smooth nbranches on LH distribution")
        test_nodes =[test_node.up for test_node in test_nodes]


