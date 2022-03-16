#!/usr/bin/env python

from tkinter.tix import TCL_TIMER_EVENTS
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
import time
import pandas as pd

from treetime import TreeTime
from treetime.utils import parse_dates
from treetime.node_interpolator import NodeInterpolator
from treetime.branch_len_interpolator import BranchLenInterpolator

from fft_tests import get_tree_events, compare, get_large_differences, get_test_node

def sort(df_list):
    import pandas as pd
    for i in range(0, len(df_list)):
        df_list[i] = pd.DataFrame(df_list[i]).set_index(1)
    for i in range(1, len(df_list)):
        df_list[i]= df_list[i].reindex(index=df_list[0].index)
    return df_list


def bland_altman_plot(data1, data2, name1, name2, title, ylim, *args, **kwargs):
    data1     = np.asarray(data1)
    data2     = np.asarray(data2)
    mean      = np.mean([data1, data2], axis=0)
    diff      = data1 - data2                   # Difference between data1 and data2
    md        = np.mean(diff)                   # Mean of the difference
    sd        = np.std(diff, axis=0)            # Standard deviation of the difference

    plt.rcParams['font.size'] = '6'
    fig = plt.scatter(mean, diff, *args, **kwargs)
    plt.axhline(md,           color='gray', linestyle='--')
    plt.axhline(md + 1.96*sd, color='gray', linestyle='--')
    plt.axhline(md - 1.96*sd, color='gray', linestyle='--')
    plt.ylabel(name1 + " - " + name2, fontsize=7)
    plt.xlabel("(" + name1 + " + " + name2 + ")/2", fontsize=7)
    plt.title(title, fontsize=10)
    plt.ylim(ylim)
    return fig

def run_first_round(tree, kwargs):
    tree._set_branch_length_mode(kwargs["branch_length_mode"])
    tree.infer_ancestral_sequences(infer_gtr=False, marginal=kwargs["marginal_sequences"])
    tree.prune_short_branches()
    tree.clock_filter(reroot='least-squares', n_iqd=1, plot=False, fixed_clock_rate=kwargs["clock_rate"])
    tree.reroot(root='least-squares', clock_rate=kwargs["clock_rate"])
    tree.infer_ancestral_sequences(**kwargs)
    tree.make_time_tree(clock_rate=kwargs["clock_rate"], time_marginal=kwargs["time_marginal"], divide=False)



if __name__ == '__main__':

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

    dates = parse_dates(base_name+'.metadata.csv')
    ##time the calculation and assess accuracy as difference from tt_numerical and other fft
    start = time.process_time()
    tt_numerical = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', use_fft=False,
                    aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)
    kwargs = {"marginal_sequences":True,
                "branch_length_mode": 'input',
                "sample_from_profile":"root",
                "reconstruct_tip_states":False,
                "max_iter":1,
                'clock_rate':clock_rate,
                'time_marginal':'assign'}
    run_first_round(tt_numerical, kwargs)
    tree_make_time = time.process_time() - start
    tree_events_numerical = get_tree_events(tt_numerical)
    print("Time to make tree using ultra fine numerical grid: " + str(tree_make_time))

    precision_fft = [5, 25, 50, 75, 100, 150, 200, 300, 400]
    branch_grid_size = [50, 75, 100, 150, 200, 300, 400]
    time_fft = []
    time_grid_fft_list = []
    divergence_time_diff_list = []

    for b in branch_grid_size:
        pre_time_grid_fft_list =[]
        for prec in precision_fft:
            start = time.process_time()
            tt_fft = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', precision_fft=prec, use_fft=True,
                        aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, precision_branch=b, debug=True)

            tt_fft.run(root="best", branch_length_mode='input', time_marginal='assign', max_iter=1)
            tree_make_time_fft = time.process_time() - start
            print("Time to make tree using fft grid of size " + str(prec) +" :"+str(tree_make_time_fft))
            tree_events_fft = get_tree_events(tt_fft)
            if len([t[1] for t in tree_events_fft if np.isnan(t[0])])>0:
                print("error in tree! a clade has no time before present!")
            time_fft.append(tree_make_time_fft)
            pre_time_grid_fft_list.append(tree_events_fft)
            output_comparison = compare(tree_events_numerical, tree_events_fft)
            groups = output_comparison[1].groupby('bad_branch')
            plt.figure()
            for name, group in groups:
                plt.plot(output_comparison[1].time, output_comparison[1].difference, marker='o', linestyle='', ms=1, label=name)
            plt.xlabel("nodes ranging from root at 0 to most recent")
            plt.ylabel("difference time_before_present coalescent branch - master branch")
            plt.show()
            plt.savefig("Differences branch grid " + str(b) + " fft grid " + str(prec)+".png")

        sorted_treelist = sort(pre_time_grid_fft_list)
        divergence_times = []
        for i in range(0, len(sorted_treelist)):
            if len([t for t in sorted_treelist[i][0] if np.isnan(t)])>0:
                print("error in sort function! a clade's time has been lost!")
            divergence_times.append(np.array(sorted_treelist[i][0]))
        time_grid_fft_list.append([divergence_times])
        divergence_time_diff = []
        for i in range(1, len(divergence_times)):
            divergence_time_diff.append(np.mean(np.abs(divergence_times[i] - divergence_times[i-1])/(precision_fft[i] - precision_fft[i-1])))
        divergence_time_diff_list.append([divergence_time_diff])

    #for each possible branch length grid size plot the divergence time differences for different fft grid sizes
    plt.figure()
    for i in range(0, len(divergence_time_diff_list)):
        epsilon = 1e-9
        plt.plot(precision_fft[1:len(time_grid_fft_list[i][0])], divergence_time_diff_list[i][0], label=branch_grid_size[i])
    plt.legend(title="precision")
    plt.axhline(y=epsilon, color='black', linestyle='--')
    plt.ylabel("average rate of change in infered divergence time")
    plt.xlabel("FFT grid size")
    plt.title("Change in Infered Divergence Times - Average Rate of Change")
    plt.savefig("Optimal_Grid_Size.png")

    time_branch = []
    pre_time_grid_fft_list =[]
    for b in branch_grid_size:
        start = time.process_time()
        tt_fft = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', precision_fft=200, use_fft=True,
                    aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, precision_branch = b, debug=True)

        run_first_round(tt_fft, kwargs)
        fig, axs = plt.subplots(1,2, sharey=True, figsize=(12,8))
        Phylo.draw(tt_numerical.tree, label_func=lambda x:"", axes=axs[0])
        axs[0].set_title("numerical time tree")
        Phylo.draw(tt_fft.tree, label_func=lambda x:"", axes=axs[1])
        axs[1].set_title("fft time tree " + str(b))
        tree_make_time_fft = time.process_time() - start
        print("Time to make tree using fft grid of size " + str(150) +" :"+str(tree_make_time_fft))
        tree_events_fft = get_tree_events(tt_fft)
        if len([t[1] for t in tree_events_fft if np.isnan(t[0])])>0:
            print("error in tree! a clade has no time before present!")
        time_branch.append(tree_make_time_fft)
        pre_time_grid_fft_list.append(tree_events_fft)
        output_comparison = compare(tree_events_numerical, tree_events_fft)
        groups = output_comparison[1].groupby('bad_branch')
        plt.figure()
        for name, group in groups:
            plt.plot(output_comparison[1].time, output_comparison[1].difference, marker='o', linestyle='', ms=1, label=name)
        plt.xlabel("nodes ranging from root at 0 to most recent")
        plt.ylabel("difference time_before_present coalescent branch - master branch")
        plt.title(str(b))
        plt.show()

        large_differences = get_large_differences(output_comparison[1])

        # for n in list(large_differences.index):
        #     test_node, test_node_fft = get_test_node([tt_numerical, tt_fft], n)
        #     print([c.name for c in test_node.clades])
        #     print([c.name for c in test_node_fft.clades])
        #     if test_node.name != tt_numerical.tree.root.name and not test_node.marginal_pos_LH.is_delta:
        #         plt.figure()
        #         plt.plot(test_node.marginal_pos_LH.x, test_node.marginal_pos_LH.prob_relative(test_node.marginal_pos_LH.x), marker='o', linestyle='', ms=2, label='numerical')
        #         plt.plot(test_node_fft.marginal_pos_LH.x, test_node_fft.marginal_pos_LH.prob_relative(test_node_fft.marginal_pos_LH.x), marker='o', linestyle='', ms=1, label='fft')
        #         plt.xlim((test_node.time_before_present-0.0001,test_node.time_before_present+0.0001))
        #         plt.title(test_node.name + "LH")
        #         plt.legend()
        #         plt.show()

    sorted_treelist = sort(pre_time_grid_fft_list)
    divergence_times = []
    for i in range(0, len(sorted_treelist)):
        divergence_times.append(np.array(sorted_treelist[i][0]))
        if len([t for t in sorted_treelist[i][0] if np.isnan(t)])>0:
            print("error in sort function! a clade's time has been lost!")
    divergence_time_diff = []
    for i in range(1, len(divergence_times)):
        divergence_time_diff.append(np.nanmean(np.abs(divergence_times[i] - divergence_times[i-1])/(branch_grid_size[i] - branch_grid_size[i-1])))

    plt.figure()
    epsilon = 1e-9
    plt.plot(branch_grid_size[1:], divergence_time_diff, label='rate of change')
    plt.plot(branch_grid_size, np.array(time_branch)*(max(divergence_time_diff)/max(time_branch)), label='time')
    plt.axhline(y=epsilon, color='black', linestyle='--')
    plt.ylabel("average rate of change in infered divergence time")
    plt.xlabel("branch grid size")
    plt.title("Change in Infered Divergence Times - Average Rate of Change")
    plt.legend()
    plt.savefig("Optimal_Branch_Grid_Size.png")

    ##after visualizing the differences in accuracy for all internal nodes I will also measure the difference in execution times by
    # running a conolution integral of a branch length interpolation distribution object and a likelihood distribution

    # def time_convolution(calc_type, grid_size):
    #     node = names_fft[key]
    #     tt = tt_fft
    #     bl = BranchLenInterpolator(node, tt.gtr, pattern_multiplicity = tt.data.multiplicity, min_width= tt.min_width,one_mutation=tt.one_mutation, branch_length_mode=tt.branch_length_mode)
    #     h = node.marginal_pos_LH
    #     if calc_type=="fft":
    #         start = time.process_time()
    #         res = NodeInterpolator.convolve_fft(h, bl, fft_grid_size=int(grid_size/50))
    #         conv_time = time.process_time() - start
    #         true_grid_size = len(res.y)
    #     if calc_type=="num":
    #         start = time.process_time()
    #         res = NodeInterpolator.convolve(h, bl, n_grid_points=grid_size, n_integral=grid_size)[0]
    #         conv_time = time.process_time() - start
    #         true_grid_size = len(res.y)
    #     return conv_time, true_grid_size


    # grid_size = [50, 75, 100, 150, 200, 300, 400, 600]
    # time_fft = np.empty((len(grid_size),2))
    # time_num = np.empty((len(grid_size),2))
    # for t in range(len(grid_size)):
    #     #print(t)
    #     time_fft[t] = time_convolution("fft", grid_size[t])
    #     time_num[t] = time_convolution("num", grid_size[t])


    # #now additionally plot these differences
    # fig = plt.figure()
    # plt.autoscale()
    # plt.plot(time_fft[:,1], time_fft[:,0], linestyle='-',marker='o', color='red', label='FFT')
    # plt.plot(time_num[:,1], time_num[:,0], linestyle='-',marker='o', color='green', label='numerical')
    # plt.xlabel("Size of Grid [convolution output size]", fontsize=8)
    # plt.ylabel("Calculation Time [sec]", fontsize=8)
    # plt.title("Calculation Speed: " + str(names_numerical_3[key].name), fontsize=10)
    # plt.legend(prop={"size":8})
    # #plt.show()
    # pp.savefig(fig)
    # plt.close()
    # pp.close()

