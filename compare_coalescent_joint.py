import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib import cm

from treetime import TreeTime
from treetime_fft import TreeTime as TreeTimeFFT
from treetime.utils import parse_dates


def get_test_nodes(tt1, tt2, pattern):
    return [n for n in tt1.tree.find_clades() if pattern in n.name][0], [n for n in tt2.tree.find_clades() if pattern in n.name][0]

def lookup_by_names(tree):
    names = {}
    for clade in tree.get_clades():
        if clade.name != tree.root.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names


if __name__ == '__main__':

    # load data and parse dates
    plt.ion()
    base_name = '../treetime_examples/data/h3n2_na/h3n2_na_20'

    dates = parse_dates(base_name+'.metadata.csv')
    tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk', use_fft=False,
                  aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)
    tt_fft = TreeTimeFFT(gtr='Jukes-Cantor', tree = base_name+'.nwk', use_fft=False,
                  aln = base_name+'.fasta', verbose = 1, dates = dates, precision=3, debug=True)

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


    terminals = [clade.name for clade in tt.tree.get_terminals()]

    node_pattern = 'NODE_0000003'
    test_node, test_node_fft = get_test_nodes(tt, tt_fft, node_pattern)

    plt.figure()
    t = np.linspace(0,4*fixed_clock_rate,1000)
    plt.plot(t, test_node.branch_length_interpolator.prob_relative(t), label='old', ls='-')
    plt.plot(t, test_node_fft.branch_length_interpolator.prob_relative(t), label='new', ls='--')
    plt.savefig("test.png")

    Tc=0.01
    tt.add_coalescent_model(Tc)
    tt_fft.add_coalescent_model(Tc)

    plt.figure()
    plt.plot(t, tt.merger_model.cost(test_node.time_before_present, t))
    plt.plot(t, tt_fft.merger_model.cost(test_node.time_before_present, t))
    plt.savefig("test.png")

    plt.figure()
    plt.plot(t, tt.merger_model.integral_merger_rate(t+0.02))
    plt.plot(t, tt_fft.merger_model.integral_merger_rate(t+0.02))
    plt.savefig("test.png")


    plt.figure()
    t = np.linspace(0,4*fixed_clock_rate,1000)
    plt.plot(t, test_node.branch_length_interpolator.prob_relative(t), label='old', ls='-')
    plt.plot(t, test_node_fft.branch_length_interpolator.prob_relative(t), label='new', ls='--')
    plt.savefig("test.png")

    def undo_merger(test_node, t):
        return test_node.branch_length_interpolator.prob(t)*np.exp(tt.merger_model.cost(test_node.time_before_present, t))

    node_pattern = 'Indiana'
    test_node, test_node_fft = get_test_nodes(tt, tt_fft, node_pattern)
    while test_node:
        if test_node.name != tt.tree.root.name:
            plt.figure()
            t = np.linspace(0,4*fixed_clock_rate,1000)
            plt.plot(t, undo_merger(test_node, t), label='old', ls='-')
            plt.plot(t, test_node_fft.branch_length_interpolator.prob(t), label='new', ls='--')
            plt.savefig("test.png")
            print(np.max(abs(undo_merger(test_node, t)- test_node_fft.branch_length_interpolator.prob(t))))
        test_node, test_node_fft = test_node.up, test_node_fft.up

    ##the greatest absolute difference is 5.551115123125783e-16

    tt.init_date_constraints(clock_rate=fixed_clock_rate)
    tt_fft.init_date_constraints(clock_rate=fixed_clock_rate)
    tt._ml_t_joint()
    tt_fft._ml_t_joint()

    for pattern in terminals:
        test_node, test_node_fft = get_test_nodes(tt, tt_fft, pattern)
        plt.figure()
        t = np.linspace(0,4*fixed_clock_rate,1000)
        plt.plot(t, undo_merger(test_node, t), label='old', ls='-')
        plt.plot(t, test_node_fft.branch_length_interpolator.prob(t), label='new', ls='--')
        plt.show()
        plt.savefig("test.png")
        print(np.max(abs(undo_merger(test_node, t)- test_node_fft.branch_length_interpolator.prob(t))))

    ##Lx is correctly set on the terminals
    def undo_C(test_node, t):
        multiplicity = len(test_node.up.clades)
        res = test_node.joint_pos_Lx.prob_relative(t)*np.exp(-tt_fft.merger_model.integral_merger_rate(t))*(tt_fft.merger_model.total_merger_rate(t)**((multiplicity-1)/(multiplicity)))
        return res/res.max()

    node_pattern = 'Indiana'
    test_node, test_node_fft = get_test_nodes(tt, tt_fft, node_pattern)
    from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('coalescent_k_joint.pdf')
    while test_node:
        print(test_node.time_before_present - test_node_fft.time_before_present)
        if test_node.name != tt.tree.root.name:
            t_node = np.linspace(test_node.time_before_present-0.005,test_node.time_before_present+0.005,1000)
            plt.figure()
            plt.title(f"joint posterior Cx {test_node.name}, mult: {len(test_node.clades)}")
            plt.plot(t_node, test_node.joint_pos_Cx.prob_relative(t_node), label='old', ls='-')
            plt.plot(t_node, test_node_fft.joint_pos_Cx.prob_relative(t_node), label='new', ls='-')
            #plt.plot(t_node, undo_C(test_node_fft, t_node), label='new_corrected', ls='--')
            plt.legend()
            plt.show()
            pp.savefig()
            plt.savefig("test.png")
            plt.close()
            t_node = np.linspace(test_node.time_before_present-0.005,test_node.time_before_present+0.01,1000)
        plt.figure()
        plt.title(f"joint posterior Lx {test_node.name}, mult: {len(test_node.clades)}")
        plt.plot(t_node, test_node.joint_pos_Lx.prob_relative(t_node), label='old', ls='-')
        plt.plot(t_node, test_node_fft.joint_pos_Lx.prob_relative(t_node), label='new', ls='-')
        if test_node.name != tt.tree.root.name:
            plt.plot(t_node, undo_C(test_node_fft, t_node), label='new_corrected', ls='--')
        plt.legend()
        plt.show()
        pp.savefig()
        plt.savefig("test.png")
        plt.close()

        test_node, test_node_fft = test_node.up, test_node_fft.up
    pp.close()