#!/usr/bin/env python
 
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
import time

from treetime_fft import TreeTime
from treetime_fft.utils import parse_dates
from treetime_fft.node_interpolator import NodeInterpolator
from treetime_fft.branch_len_interpolator import BranchLenInterpolator


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

def lookup_by_names(tree):
    names = {}
    for clade in tree.get_nonterminals():
        if clade.name != tree.root.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names


if __name__ == '__main__':

    sc2 = False
    print("running marginal distribution")
    # load data and parse dates
    plt.ion()
    if sc2:
        tree = '/scicore/home/neher/parann00/SC2_data/tree_raw.nwk'
        #aln = '../../treetime_examples/data/SC2_data/sequences.fasta'
        aln = None
        seq_len=29903
        dates = parse_dates('/scicore/home/neher/parann00/SC2_data/metadata.tsv')
        root = 'Wuhan/Hu-1/2019'
    else:
        base_name = './treetime_examples/data/h3n2_na/h3n2_na_20'
        tree = base_name+'.nwk'
        aln = base_name+'.fasta'
        seq_len = None
        dates = parse_dates(base_name+'.metadata.csv')
        root = None

    ##I compute the tree using three different calculations of the marginal distribution, the old numerical approach 
    # (either using the default or an ultra fine grid spacing) and the new approach using the FFT. 
    
    #start = time.process_time()
    #tt_numerical_1 = TreeTime(gtr='Jukes-Cantor', tree = tree,
    #              aln = aln, verbose = 1, dates = dates, root=root, seq_len =seq_len, debug=True, precision=1, use_fft=False)

    #tt_numerical_1.run(root="best", branch_length_mode='marginal', time_marginal=True, max_iter=2)
    #tree_make_time_1 = time.process_time() - start
    #print("Time to make tree using default numerical grid: " + str(tree_make_time_1))
    start = time.process_time()
    tt_numerical_3 = TreeTime(gtr='Jukes-Cantor', tree = tree,
                  aln = aln, verbose = 1, dates = dates, root=root, seq_len =seq_len, debug=True, use_fft=False)

    tt_numerical_3.run(root="best", branch_length_mode='marginal', time_marginal=True, max_iter=2)
    tree_make_time_3 = time.process_time() - start
    print("Time to make tree using ultra fine numerical grid: " + str(tree_make_time_3))

    start = time.process_time()
    tt_fft = TreeTime(gtr='Jukes-Cantor', tree = tree,
                  aln = aln, verbose = 1, dates = dates, root=root, seq_len =seq_len, debug=True, use_fft=True)

    tt_fft.run(root="best", branch_length_mode='marginal', time_marginal=True, max_iter=2)
    tree_make_time_fft = time.process_time() - start
    print("Time to make tree using default fft grid: " + str(tree_make_time_fft))

    start = time.process_time()
    tt_fft_400 = TreeTime(gtr='Jukes-Cantor', tree = tree,
                  aln = aln, verbose = 1, dates = dates, root=root, seq_len =seq_len, debug=True, precision_fft=400, use_fft=True)

    tt_fft_400.run(root="best", branch_length_mode='marginal', time_marginal=True, max_iter=2)
    tree_make_time_fft4 = time.process_time() - start
    print("Time to make tree using 2x default fft grid: " + str(tree_make_time_fft4))
    
    #names_numerical_1 = lookup_by_names(tt_numerical_1.tree)
    names_numerical_3 = lookup_by_names(tt_numerical_3.tree)
    names_fft = lookup_by_names(tt_fft.tree)
    names_fft_4 = lookup_by_names(tt_fft_400.tree)

    def get_x_y(node, treetime_obj):
        x = node.marginal_pos_LH.x/treetime_obj.date2dist.clock_rate
        y = np.exp(-node.marginal_pos_LH.y+min(node.marginal_pos_LH.y))
        return x, y
    
    pp = PdfPages('/scicore/home/neher/parann00/treetime/test/pictures/FFT_Marginal_Likelihood_Plots.pdf')
    from matplotlib import rcParams
    rcParams.update({'figure.autolayout': True})

    firstPage = plt.figure()
    firstPage.clf()
    txt = 'Time for ultra fine numerical grid: ' + str(tree_make_time_3) + '\n Time for default fft grid: ' +str(tree_make_time_fft) + '\n Time for 2x default fft grid: ' + str(tree_make_time_fft4)
    firstPage.text(0.5,0.5,txt, transform=firstPage.transFigure, size=10, ha="center")
    pp.savefig()
    plt.show()
    plt.close()

    node_number = 0
    for key in names_numerical_3:
        node_number +=1 
        if (node_number >20):
            break
        ##calculate the x and y values for plotting
        #x_num_1, y_num_1 = get_x_y(names_numerical_1[key], tt_numerical_1)
        #print("numeric default, N: " + str(len(x_num_1)) + " N^2: " + str(len(x_num_1)**2))
        x_num_3, y_num_3 = get_x_y(names_numerical_3[key], tt_numerical_3)
        num_size = "numeric default, N: " + str(len(x_num_3)) + " N^2: " + str(round(len(x_num_3)**2))
        x_fft, y_fft = get_x_y(names_fft[key], tt_fft)
        fft_size = "FFT, N: " + str(len(x_fft)) + " 4Nlog(2N): " + str(round(4*len(x_fft)*np.log2(2*len(x_fft))))
        x_fft_4, y_fft_4 = get_x_y(names_fft_4[key], tt_fft_400)
        fft4_size = "FFT double, N: " + str(len(x_fft_4)) + " 4Nlog(2N): " + str(round(4*len(x_fft_4)*np.log2(2*len(x_fft_4))))
        ##get the effective support (calculated by fft version), using for making xlim the same for all plots
        effective_support = names_fft_4[key].marginal_pos_LH.effective_support/tt_fft_400.date2dist.clock_rate
        # if effective_support[0]< max(x_num_1[0], x_num_3[0],x_fft[0], x_fft_4[0]):
        #     print("effective support not fully covered, old start point: "+ str(effective_support[0]))
        #     effective_support[0] = max(x_num_1[0], x_num_3[0],x_fft[0], x_fft_4[0])
        #     print("new start point: "+ str(effective_support[0]))
        if effective_support[0]< max(x_num_3[0],x_fft[0], x_fft_4[0]):
            print("effective support not fully covered, old start point: "+ str(effective_support[0]))
            effective_support[0] = max(x_num_3[0],x_fft[0], x_fft_4[0])
            print("new start point: "+ str(effective_support[0]))
        eff_support_uniform = np.linspace(effective_support[0], effective_support[1], num=1000, endpoint=True)

        firstPage = plt.figure()
        firstPage.clf()
        txt = num_size + '\n' + fft_size + '\n' + fft4_size
        firstPage.text(0.5,0.5,txt, transform=firstPage.transFigure, size=10, ha="center")
        pp.savefig()
        plt.close()

        ##plot the three distributions
        fig = plt.figure()
        #fig = plt.figure()
        #plt.autoscale()
        #plt.plot(x_num_1, y_num_1, 'o', markersize=1, color='red', label='numerical (default))')
        plt.plot(x_num_3, y_num_3, 'o', markersize=1, color='green', label='numerical (ultra fine grid))')
        plt.plot(x_fft, y_fft, 'o', markersize=1, color='blue', label='FFT (default FWHM grid)')
        plt.plot(x_fft_4, y_fft_4, 'o', markersize=1, color='black', label='FFT (double default FWHM grid)')
        plt.xlim(effective_support)
        plt.xlabel("Time (effective support)", fontsize=8)
        plt.ylabel("P", fontsize=8)
        plt.title("Marginal Likelihood (Time): " + str(names_numerical_3[key].name), fontsize=10)
        plt.legend(prop={"size":8})
        #plt.show()
        pp.savefig(fig)
        plt.close()

        #numer_1 = interp1d(x_num_1, y_num_1)
        numer_3 = interp1d(x_num_3, y_num_3)
        fft = interp1d(x_fft, y_fft)
        fft_4 = interp1d(x_fft_4, y_fft_4)

        # ##perform a Bland-Altman plot of the distributions (comparing the difference in accuracy of the default FFT with two different numerical grids)
        
        
        # ylim_min = min(min(numer_1(eff_support_uniform)-fft(eff_support_uniform)), min(numer_3(eff_support_uniform)-fft(eff_support_uniform)))
        # ylim_max = max(max(numer_1(eff_support_uniform)-fft(eff_support_uniform)), max(numer_3(eff_support_uniform)-fft(eff_support_uniform)))
        # fig = plt.figure()
        # plt.rcParams['font.size'] = '6'
        # plt.autoscale()
        # plt.subplot(1,2,1)
        # fig.tight_layout()
        # bland_altman_plot(numer_1(eff_support_uniform), fft(eff_support_uniform), "numerical calculation", "FFT calculation", "default grid: " +str(names_numerical_1[key].name), ylim=(ylim_min, ylim_max))
        # sp = plt.subplot(1,2,2)
        # sp.yaxis.tick_right()
        # sp.yaxis.set_label_position("right")
        # bland_altman_plot(numer_3(eff_support_uniform), fft(eff_support_uniform), "numerical calculation", "FFT calculation", "ultra fine grid: " +str(names_numerical_1[key].name), ylim=(ylim_min, ylim_max))
        # pp.savefig(fig, bb_inches="tight")

        # ##perform a Bland-Altman plot of the distributions (comparing the difference in accuracy of two FFT grids with default numerical grids)
        # ylim_min = min(min(numer_1(eff_support_uniform)-fft(eff_support_uniform)), min(numer_1(eff_support_uniform)-fft_4(eff_support_uniform)))
        # ylim_max = max(max(numer_1(eff_support_uniform)-fft(eff_support_uniform)), max(numer_1(eff_support_uniform)-fft_4(eff_support_uniform)))
        # fig = plt.figure()
        # plt.rcParams['font.size'] = '6'
        # plt.autoscale()
        # plt.subplot(1,2,1)
        # fig.tight_layout()
        # bland_altman_plot(numer_1(eff_support_uniform), fft(eff_support_uniform), "numerical calculation", "FFT calculation", "default FWHM grid: " +str(names_numerical_1[key].name), ylim=(ylim_min, ylim_max))
        # sp = plt.subplot(1,2,2)
        # sp.yaxis.tick_right()
        # sp.yaxis.set_label_position("right")
        # bland_altman_plot(numer_1(eff_support_uniform), fft_4(eff_support_uniform), "numerical calculation", "FFT calculation", "2x default FWHM grid: " +str(names_numerical_1[key].name), ylim=(ylim_min, ylim_max))
        # pp.savefig(fig, bb_inches="tight")
        ##notably the difference between the calculations increases at a higher rate when the grid size of the FFT grid is doubled than when the grid size of the numerical approach is doubled
        #therefore it would make sense to check that doubling the grid size in both approaches makes the difference between them smaller (i.e. that the two solutions actually converge)
        
        
        ylim_min = min(min(numer_3(eff_support_uniform)-fft(eff_support_uniform)), min(numer_3(eff_support_uniform)-fft_4(eff_support_uniform)))
        ylim_max = max(max(numer_3(eff_support_uniform)-fft(eff_support_uniform)), max(numer_3(eff_support_uniform)-fft_4(eff_support_uniform)))
        fig = plt.figure()
        plt.rcParams['font.size'] = '6'
        plt.autoscale()
        plt.subplot(1,2,1)
        fig.tight_layout()
        bland_altman_plot(numer_3(eff_support_uniform), fft(eff_support_uniform), "numerical calculation", "FFT calculation", "default grids: " +str(names_numerical_3[key].name), ylim=(ylim_min, ylim_max))
        sp = plt.subplot(1,2,2)
        sp.yaxis.tick_right()
        sp.yaxis.set_label_position("right")
        bland_altman_plot(numer_3(eff_support_uniform), fft_4(eff_support_uniform), "numerical calculation", "FFT calculation", "2x default grids: " +str(names_numerical_3[key].name), ylim=(ylim_min, ylim_max))
        pp.savefig(fig, bb_inches="tight")
        plt.close()
        ##there is still a notable difference in the output of these methods and convergence is not as clear. 
        ##It would make sense to look at functions where we know the analytical result of a convolution to be able to properly assess the accuracy of our methods
        
        eff_support_uniform_right = np.linspace(effective_support[0], effective_support[1]+1, num=1000, endpoint=True)
        ##plot normalized difference along time axis
        epsilon = 1e-11
        plt.rcParams['font.size'] = '6'
        fig = plt.figure()
        plt.autoscale()
        ##get the maximum and minimum values of both normalized differences for ylim so that plots are comparable
        y1 = (numer_3(eff_support_uniform_right) - fft(eff_support_uniform_right))/ (numer_3(eff_support_uniform_right)+epsilon)
        y3 = (numer_3(eff_support_uniform_right) - fft_4(eff_support_uniform_right))/ (numer_3(eff_support_uniform_right)+epsilon)
        ylim_min = min(min(y1), min(y3))
        ylim_max = max(max(y1), max(y3))
        plt.subplot(1,2,1)
        fig.tight_layout()
        plt.scatter(eff_support_uniform_right, y1)
        plt.xlabel("Time (effective support)", fontsize=7)
        plt.ylabel("(numerical calculation - FFT calculation)/numerical calculation", fontsize=7)
        plt.ylim((ylim_min, ylim_max))
        plt.title("default grid: " + str(names_numerical_3[key].name), fontsize=10)
        sp = plt.subplot(1,2,2)
        sp.yaxis.tick_right()
        sp.yaxis.set_label_position("right")
        plt.scatter(eff_support_uniform_right, y3)
        plt.xlabel("Time (effective support + 1T)", fontsize=7)
        plt.ylabel("(numerical calculation - FFT calculation)/ numerical calculation", fontsize=7)
        plt.ylim((ylim_min, ylim_max))
        plt.title("2x default grid: " + str(names_numerical_3[key].name), fontsize=10)
        pp.savefig(fig)
        plt.close()
        
        #eff_support_uniform_left = np.linspace(effective_support[0]-1, effective_support[1], num=1000, endpoint=True)
        #eff_support_uniform_right = np.linspace(effective_support[0], effective_support[1]+1, num=1000, endpoint=True)
        #eff_support_uniform_left_only = np.linspace(effective_support[0]-1, effective_support[0], num=100, endpoint=False)
        eff_support_uniform_right_only = np.linspace(effective_support[1], effective_support[1]+1, num=100, endpoint=False)
        #supp = [eff_support_uniform_left_only, eff_support_uniform_left, eff_support_uniform, eff_support_uniform_right, eff_support_uniform_right_only]
        #names_supp = ['effective support, 1T earlier', 'effective support, shifted 1T earlier', 'effective support', 'effective support, shifted 1T later', 'effective support, 1T later']
        supp = [ eff_support_uniform, eff_support_uniform_right_only]
        names_supp = ['effective support', '1T after effective support']
        ##plot the logarithm of the three distributions
        fig = plt.figure()
        for i in range(len(supp)):
            plt.subplot(1,2,i+1)
            eff_support_uniform_new = supp[i]
            plt.autoscale()
            #plt.plot(eff_support_uniform_new, -names_numerical_1[key].marginal_pos_LH(eff_support_uniform_new*tt_fft.date2dist.clock_rate)+min(names_numerical_1[key].marginal_pos_LH.y), 'o', markersize=1, color='red', label='numerical (default))')
            plt.plot(eff_support_uniform_new, -names_numerical_3[key].marginal_pos_LH(eff_support_uniform_new*tt_fft.date2dist.clock_rate)+min(names_numerical_3[key].marginal_pos_LH.y), 'o', markersize=1, color='green', label='numerical (ultra fine grid))')
            plt.plot(eff_support_uniform_new, -names_fft[key].marginal_pos_LH(eff_support_uniform_new*tt_fft.date2dist.clock_rate)+min(names_fft[key].marginal_pos_LH.y), 'o', markersize=1, color='blue', label='FFT')
            plt.xlabel("Time (" + names_supp[i] + ")", fontsize=8)
            plt.ylabel("log(P)", fontsize=8)
            plt.title("Log(Marginal P) (Time): " + str(names_numerical_3[key].name), fontsize=10)
            plt.legend(prop={"size":8})
        #plt.show()
        pp.savefig(fig)
        plt.close()


    ##after visualizing the differences in accuracy for all internal nodes I will also measure the difference in execution times by 
    # running a conolution integral of a branch length interpolation distribution object and a likelihood distribution
    
    def time_convolution(calc_type, grid_size):
        node = names_fft[key]
        tt = tt_fft
        bl = BranchLenInterpolator(node, tt.gtr, pattern_multiplicity = tt.data.multiplicity, min_width= tt.min_width,one_mutation=tt.one_mutation, branch_length_mode=tt.branch_length_mode)
        h = node.marginal_pos_LH
        if calc_type=="fft":
            start = time.process_time()
            res = NodeInterpolator.convolve_fft(h, bl, fft_grid_size=int(grid_size/50))
            conv_time = time.process_time() - start
            true_grid_size = len(res.y)
        if calc_type=="num":
            start = time.process_time()
            res = NodeInterpolator.convolve(h, bl, n_grid_points=grid_size, n_integral=grid_size)[0]
            conv_time = time.process_time() - start
            true_grid_size = len(res.y)
        return conv_time, true_grid_size

    
    grid_size = [50, 75, 100, 150, 200, 300, 400, 600]
    time_fft = np.empty((len(grid_size),2))
    time_num = np.empty((len(grid_size),2))
    for t in range(len(grid_size)):
        #print(t)
        time_fft[t] = time_convolution("fft", grid_size[t])
        time_num[t] = time_convolution("num", grid_size[t])


    #now additionally plot these differences
    fig = plt.figure()
    plt.autoscale()
    plt.plot(time_fft[:,1], time_fft[:,0], linestyle='-',marker='o', color='red', label='FFT')
    plt.plot(time_num[:,1], time_num[:,0], linestyle='-',marker='o', color='green', label='numerical')
    plt.xlabel("Size of Grid [convolution output size]", fontsize=8)
    plt.ylabel("Calculation Time [sec]", fontsize=8)
    plt.title("Calculation Speed: " + str(names_numerical_3[key].name), fontsize=10)
    plt.legend(prop={"size":8})
    #plt.show()
    pp.savefig(fig)
    plt.close()
    pp.close()


