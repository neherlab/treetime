#!/usr/bin/env python

## code for determining the optimal grid spacing for the FFT version, using time and memory usage
## as well as accuracy estimations based on convergence and positional likelihood of tree
## needs an output results folder

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import time
import pandas as pd
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

from treetime import TreeTime
from treetime.utils import parse_dates
from resource import getrusage, RUSAGE_SELF
from fft_tests import get_tree_events, compare, plot_differences

def sort(df_list):
    for i in range(0, len(df_list)):
        df_list[i] = pd.DataFrame(df_list[i]).set_index(1)
    for i in range(1, len(df_list)):
        df_list[i]= df_list[i].reindex(index=df_list[0].index)
    return df_list


def run_first_round(tree, kwargs):
    tree._set_branch_length_mode(kwargs["branch_length_mode"])
    tree.infer_ancestral_sequences(infer_gtr=False, marginal=kwargs["marginal_sequences"])
    tree.prune_short_branches()
    tree.clock_filter(reroot='least-squares', n_iqd=1, plot=False, fixed_clock_rate=kwargs["clock_rate"])
    tree.reroot(root='least-squares', clock_rate=kwargs["clock_rate"])
    tree.infer_ancestral_sequences(**kwargs)
    tree.make_time_tree(clock_rate=kwargs["clock_rate"], time_marginal=kwargs["time_marginal"], divide=False)
    peak_memory = getrusage(RUSAGE_SELF).ru_maxrss
    return peak_memory

def plot_accuracy_memory_time(branch_grid_size, accuracy, time_branch, memory_branch, accuracy_label, title):
    plt.figure()

    host = host_subplot(111, axes_class=AA.Axes)
    plt.subplots_adjust(right=0.75)

    par1 = host.twinx()
    par2 = host.twinx()

    offset = 60
    new_fixed_axis = par2.get_grid_helper().new_fixed_axis
    par2.axis["right"] = new_fixed_axis(loc="right", axes=par2, offset=(offset, 0))
    par2.axis["right"].toggle(all=True)
    par1.axis["right"].toggle(all=True)

    host.set_xlabel("Branch Grid Size")
    host.set_ylabel(accuracy_label)
    par1.set_ylabel("Time [sec]")
    par2.set_ylabel("Memory Consumption [MB]")

    if len(accuracy)==(len(branch_grid_size)-1):
        p1, = host.plot(branch_grid_size[1:], accuracy, label=title)
        p2, = par1.plot(branch_grid_size, np.array(time_branch), label='time')
        p3, = par2.plot(branch_grid_size, np.array(memory_branch), label='memory consumption')
    else:
        p1, = host.plot(branch_grid_size, accuracy, label=title)
        p2, = par1.plot(branch_grid_size, np.array(time_branch), label='time')
        p3, = par2.plot(branch_grid_size, np.array(memory_branch), label='memory consumption')

    host.legend()

    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    par2.axis["right"].label.set_color(p3.get_color())

    plt.draw()
    plt.show()
    plt.title(title + " - FFT Grid Size=" + str(prec_fft))
    plt.savefig("./results/Memory_Time_" + str(title)+ ".png")


if __name__ == '__main__':

    ##model parameters for testing
    # choose if should be tested on ebola or h3n2_na dataset
    ebola=True

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
    memory_numerical = run_first_round(tt_numerical, kwargs)
    print(memory_numerical)
    tree_make_time = time.process_time() - start
    tree_events_numerical = get_tree_events(tt_numerical)
    print("Time to make tree using ultra fine numerical grid: " + str(tree_make_time))

    numerical_LH= tt_numerical.tree.positional_LH
    print("likelihood of tree under numerical:" + str(numerical_LH))

    precision_fft = [5, 25, 50, 75, 100, 150, 200, 300, 400]
    branch_grid_size = [50, 75, 100, 150, 200, 300, 400]

    def find_optimal_branch_fft_gridsizes(precision_fft, branch_grid_size):
        '''
        Given an input list of FFT and branch grid sizes perform the first iteration of the run function for each combination,
        determine the inferred divergence times for each new timetree and plot the average rate of change in inferred divergence time
        for each branch length when changing the fft grid size
        '''
        time_fft = [] #time needed to make fft
        divergence_times_fft_list = []
        divergence_time_diff_list = []

        pp = PdfPages('./results/DifferencesComparison.pdf') #output the difference plots to visually assess if there is an issue

        for b in branch_grid_size:
            pre_divergence_times_fft_list =[]
            pre_time_fft = []
            for prec in precision_fft:
                start = time.process_time()
                tt_fft = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', precision_fft=prec, use_fft=True,
                            aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, precision_branch=b, debug=True)

                run_first_round(tt_fft, kwargs)
                tree_make_time_fft = time.process_time() - start
                print("Time to make tree using fft grid of size " + str(prec) +" :"+str(tree_make_time_fft))
                tree_events_fft = get_tree_events(tt_fft)
                if len([t[1] for t in tree_events_fft if np.isnan(t[0])])>0:
                    print("error in tree! a clade has no time before present!")
                pre_time_fft.append(tree_make_time_fft)
                pre_divergence_times_fft_list.append(tree_events_fft)
                output_comparison = compare(tree_events_numerical, tree_events_fft)
                fig = plot_differences(output_comparison[1], title="Differences branch grid " + str(b) + " fft grid " + str(prec))
                pp.savefig(fig)
            time_fft.append([pre_time_fft])
            sorted_treelist = sort(pre_divergence_times_fft_list)
            divergence_times = []
            for i in range(0, len(sorted_treelist)):
                if len([t for t in sorted_treelist[i][0] if np.isnan(t)])>0:
                    print("error in sort function! a clade's time has been lost!")
                divergence_times.append(np.array(sorted_treelist[i][0]))
            divergence_times_fft_list.append([divergence_times])
            divergence_time_diff = []
            for i in range(1, len(divergence_times)):
                divergence_time_diff.append(np.mean(np.abs(divergence_times[i] - divergence_times[i-1])/(precision_fft[i] - precision_fft[i-1])))
            divergence_time_diff_list.append([divergence_time_diff])

        pp.close()

        #for each possible branch length grid size plot the divergence time differences for different fft grid sizes
        plt.figure()
        for i in range(0, len(divergence_time_diff_list)):
            epsilon = 1e-9
            plt.plot(precision_fft[1:len(divergence_times_fft_list[i][0])], divergence_time_diff_list[i][0], 'o-', label=branch_grid_size[i])
        plt.legend(title="branch grid size")
        plt.axhline(y=epsilon, color='black', linestyle='--')
        plt.ylabel("log(average rate of change)")
        plt.xlabel("FFT grid size")
        plt.yscale('log')
        plt.title("Change in Infered Divergence Times - Average Rate of Change")
        plt.savefig("./results/Optimal_Grid_Size.png")

        return time_fft, divergence_times_fft_list, divergence_time_diff_list
    time_fft, divergence_times_fft_list, divergence_time_diff_list = find_optimal_branch_fft_gridsizes(precision_fft, branch_grid_size)

    prec_fft = 175
    def find_optimal_branch_gridsizes(prec_fft, branch_grid_size):
        """
        For a given (optimal) FFT grid size plot the change of accuracy (either by the average rate of change or the LH)
        and the memory and time consumption at each branch grid size
        """
        time_branch = []
        memory_branch = []
        fft_LH = []

        pre_divergence_times_fft_list =[]
        for b in branch_grid_size:
            start = time.process_time()
            tt_fft = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', precision_fft=prec_fft, use_fft=True,
                        aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, precision_branch = b, debug=True)

            memory_fft = run_first_round(tt_fft, kwargs)
            tree_make_time_fft = time.process_time() - start
            print("Time to make tree using fft grid of size " + str(150) +" :"+str(tree_make_time_fft))
            tree_events_fft = get_tree_events(tt_fft)
            if len([t[1] for t in tree_events_fft if np.isnan(t[0])])>0:
                print("error in tree! a clade has no time before present!")
            time_branch.append(tree_make_time_fft)
            memory_branch.append(memory_fft)
            pre_divergence_times_fft_list.append(tree_events_fft)
            output_comparison = compare(tree_events_numerical, tree_events_fft)
            plot_differences(output_comparison[1], title="Differences branch grid " + str(b) + " fft grid " + str(prec_fft))
            fft_LH.append(tt_fft.tree.positional_LH)
            print("likelihood of tree under fft:" + str(tt_fft.tree.positional_LH))


        sorted_treelist = sort(pre_divergence_times_fft_list)
        divergence_times = []
        for i in range(0, len(sorted_treelist)):
            divergence_times.append(np.array(sorted_treelist[i][0]))
            if len([t for t in sorted_treelist[i][0] if np.isnan(t)])>0:
                print("error in sort function! a clade's time has been lost!")
        divergence_time_diff = []
        for i in range(1, len(divergence_times)):
            divergence_time_diff.append(np.nanmean(np.abs(divergence_times[i] - divergence_times[i-1])/(branch_grid_size[i] - branch_grid_size[i-1])))

        memory_branch = np.array(memory_branch)/(10**6)
        accuracy_rate = np.array(divergence_time_diff)/1e-8
        plot_accuracy_memory_time(branch_grid_size, accuracy_rate, time_branch, memory_branch, "Average Rate of Change Inferred Time [1e-8]", "average_rate_of_change")
        plot_accuracy_memory_time(branch_grid_size, np.array(fft_LH), time_branch, memory_branch, "tree log(LH) ", "log(LH)")

        return time_branch, memory_branch, divergence_time_diff, fft_LH

    time_branch, memory_branch, divergence_time_diff, fft_LH = find_optimal_branch_gridsizes(prec_fft, branch_grid_size)

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

