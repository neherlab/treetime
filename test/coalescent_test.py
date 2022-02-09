#!/usr/bin/env python
 
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
import time
from coalescent_test_helper_functions import initialize_basetree, add_coalescent, multiply

from treetime import TreeTime
from treetime.merger_models import Coalescent
from treetime.utils import parse_dates


def lookup_by_names(tree):
    names = {}
    for clade in tree.get_nonterminals():
        if clade.name != tree.root.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

def lookup_by_names_terminals(tree):
    names = {}
    for clade in tree.get_terminals():
        if clade.name != tree.root.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

def eq(self1, other):
    """"
    Check if the objects self1 and other are equal, i.e. their objects are the same
    """
    same = True
    var_dict = self1.__dict__
    #print(var_dict)
    l = ['_merger_cost', '_func', 'gtr', 'node', 'logger', 'profile_map', 'up', 'tt', 'clades']
    for key in var_dict:
        if key not in l:
            try:
                if not np.all(getattr(self1, key)==getattr(other, key)):
                    print(key)
                    same = False
            except:
                print("the key: " + str(key) + "does not work")
    return same

##define variables
base_name = '../../treetime_examples/data/h3n2_na/h3n2_na_20'
tree = base_name+'.nwk'
aln = base_name+'.fasta'
seq_len = None
dates = parse_dates(base_name+'.metadata.csv')
root = None
Tc=0.01
fixed_clock_rate=0.002

##initialize trees
tt_numerical_old = TreeTime(gtr='Jukes-Cantor', tree = tree,
                aln = aln, verbose = 1, dates = dates, root=root, seq_len =seq_len, debug=True, use_fft=False)

initialize_basetree(tt_numerical_old, root="best", branch_length_mode='marginal', time_marginal=True, max_iter=1, Tc=Tc, fixed_clock_rate=fixed_clock_rate)

from treetime_fft import TreeTime as TreeTimeFFT
tt_numerical_new = TreeTimeFFT(gtr='Jukes-Cantor', tree = tree,
                aln = aln, verbose = 1, dates = dates, root=root, seq_len =seq_len, debug=True, use_fft=False)

initialize_basetree(tt_numerical_new, root="best", branch_length_mode='marginal', time_marginal=True, max_iter=1, Tc=Tc, fixed_clock_rate=fixed_clock_rate)

##get names to facilitate looking at specific nodes in the tree
names_numerical_old = lookup_by_names(tt_numerical_old.tree)
names_numerical_new = lookup_by_names(tt_numerical_new.tree)

terminals_numerical_old = lookup_by_names_terminals(tt_numerical_old.tree)
terminals_numerical_new = lookup_by_names_terminals(tt_numerical_new.tree)


terminals = [key for key in terminals_numerical_old]
#check nodes for differences
node_number = 0
for key in terminals_numerical_old:
    if not eq(terminals_numerical_old[key], terminals_numerical_new[key]):
        print("there is an issue 1")

for key in names_numerical_old:
    if not eq(names_numerical_old[key], names_numerical_new[key]):
        print("there is an issue 1")

##if issues are flagged they can be looked at more closely using the following code
# terminals = [key for key in terminals_numerical_old]
# total_issues = 0
# for key in terminals_numerical_old:
#     if not eq(terminals_numerical_old[key].joint_pos_Cx, terminals_numerical_new[key].joint_pos_Cx):
#         print("there is an issue Cx")
#         total_issues += 1
#     if not eq(terminals_numerical_old[key].joint_pos_Lx, terminals_numerical_new[key].joint_pos_Lx):
#         print("there is an issue with Lx")
#     if not eq(terminals_numerical_old[key].branch_length_interpolator, terminals_numerical_new[key].branch_length_interpolator):
#         print("there is an issue 1")
#         total_issues += 1
#     if not eq(terminals_numerical_old[key].branch_length_interpolator.gtr, terminals_numerical_new[key].branch_length_interpolator.gtr):
#         print("there is an issue 2")
#     if not eq(terminals_numerical_old[key].branch_length_interpolator.node, terminals_numerical_new[key].branch_length_interpolator.node):
#         print("there is an issue 3")
#         if not terminals_numerical_old[key].branch_length == terminals_numerical_new[key].branch_length:
#             print("there is an issue branch length")
#             total_issues += 1
#             print(terminals_numerical_old[key].branch_length - terminals_numerical_new[key].branch_length)
#         if not terminals_numerical_old[key].branch_length_interpolator.node.clock_length == terminals_numerical_new[key].branch_length_interpolator.node.clock_length:
#             print("there is an issue clock length")
#             total_issues += 1
#             print(terminals_numerical_old[key].branch_length_interpolator.node.clock_length -terminals_numerical_new[key].branch_length_interpolator.node.clock_length)
#         if not eq(terminals_numerical_old[key].branch_length_interpolator.node.joint_pos_Lx, terminals_numerical_new[key].branch_length_interpolator.node.joint_pos_Lx):
#             print("there is an issue Lx")
#             total_issues += 1
#         if not eq(terminals_numerical_old[key].branch_length_interpolator.node.joint_pos_Cx, terminals_numerical_new[key].branch_length_interpolator.node.joint_pos_Cx):
#             print("there is an issue Cx")
#             total_issues += 1
# print("total number of issues:" +str(total_issues))

# total_issues = 0
# for key in names_numerical_old:
#     if not eq(names_numerical_old[key].joint_pos_Cx, names_numerical_new[key].joint_pos_Cx):
#         print("there is an issue Cx")
#         total_issues += 1
#     if not eq(names_numerical_old[key].joint_pos_Lx, names_numerical_new[key].joint_pos_Lx):
#         print("there is an issue with Lx")
#         total_issues += 1
#     if not eq(names_numerical_old[key].branch_length_interpolator, names_numerical_new[key].branch_length_interpolator):
#         print("there is an issue 1")
#         total_issues += 1
#     if not eq(names_numerical_old[key].branch_length_interpolator.gtr, names_numerical_new[key].branch_length_interpolator.gtr):
#         print("there is an issue 2")
#     if not eq(names_numerical_old[key].branch_length_interpolator.node, names_numerical_new[key].branch_length_interpolator.node):
#         print("there is an issue 3")
#         if not names_numerical_old[key].branch_length == names_numerical_new[key].branch_length:
#             print("there is an issue branch length")
#             total_issues += 1
#             print(names_numerical_old[key].branch_length - names_numerical_new[key].branch_length)
#         if not names_numerical_old[key].branch_length_interpolator.node.clock_length == names_numerical_new[key].branch_length_interpolator.node.clock_length:
#             print("there is an issue clock length")
#             total_issues += 1
#             print(names_numerical_old[key].branch_length_interpolator.node.clock_length -names_numerical_new[key].branch_length_interpolator.node.clock_length)
#         if not eq(names_numerical_old[key].branch_length_interpolator.node.joint_pos_Lx, names_numerical_new[key].branch_length_interpolator.node.joint_pos_Lx):
#             print("there is an issue Lx")
#             total_issues += 1
#         if not eq(names_numerical_old[key].branch_length_interpolator.node.joint_pos_Cx, names_numerical_new[key].branch_length_interpolator.node.joint_pos_Cx):
#             print("there is an issue Cx")
#             total_issues += 1
# print("total number of issues:" +str(total_issues))

##if coalescent modelnot yet added can add seperately 
# add_coalescent(tt_numerical_old, False, root="best", branch_length_mode='marginal', time_marginal=True, max_iter=1, Tc=Tc)
# add_coalescent(tt_numerical_new, True, root="best", branch_length_mode='marginal', time_marginal=True, max_iter=1, Tc=Tc)
# tt_numerical_old.merger_model.attach_to_tree()

##check that the initialized merger_models are the same
def check_merger_model():
    tt_numerical_new_tree_events = np.array(sorted([(n.time_before_present, len(n.clades)-1)
                                    for n in tt_numerical_new.tree.find_clades() if not n.bad_branch],
                                    key=lambda x:-x[0]))

    tt_numerical_old_tree_events = np.array(sorted([(n.time_before_present, len(n.clades)-1)
                                    for n in tt_numerical_old.tree.find_clades() if not n.bad_branch],
                                    key=lambda x:-x[0]))

    print(np.all(tt_numerical_new_tree_events==tt_numerical_old_tree_events))

    eq(tt_numerical_old.merger_model.Tc, tt_numerical_new.merger_model.Tc)
    print(np.isnan(tt_numerical_new.merger_model.Tc._fill_value_orig) and np.isnan(tt_numerical_old.merger_model.Tc._fill_value_orig))
    print(np.isnan(tt_numerical_new.merger_model.Tc._fill_value_below) and np.isnan(tt_numerical_old.merger_model.Tc._fill_value_below))
    print(np.isnan(tt_numerical_new.merger_model.Tc._fill_value_above) and np.isnan(tt_numerical_old.merger_model.Tc._fill_value_above))
    eq(tt_numerical_old.merger_model.integral_merger_rate, tt_numerical_new.merger_model.integral_merger_rate)
    eq(tt_numerical_old.merger_model.nbranches, tt_numerical_new.merger_model.nbranches)
    print(np.isnan(tt_numerical_new.merger_model.integral_merger_rate._fill_value_orig) and np.isnan(tt_numerical_old.merger_model.integral_merger_rate._fill_value_orig))
    print(np.isnan(tt_numerical_new.merger_model.integral_merger_rate._fill_value_below) and np.isnan(tt_numerical_old.merger_model.integral_merger_rate._fill_value_below))
    print(np.isnan(tt_numerical_new.merger_model.integral_merger_rate._fill_value_above) and np.isnan(tt_numerical_old.merger_model.integral_merger_rate._fill_value_above))

check_merger_model()


##after the merger_cost has been added to the nodes the branch length interpolation objects change, check for errors to get approximate size of numerical error
def get_new_from_old(L, node, treetime_obj, Lx=True):
    from treetime.distribution import Distribution
    multiplicity = len(node.up.clades)
    tnode = node.time_before_present
    x = L.x
    negative_merger_cost = Distribution(x, -treetime_obj.merger_model.cost(tnode, x, multiplicity), is_log=True)
    L1 = multiply(x, [L, negative_merger_cost])
    L = Distribution(x, L.__call__(x)-treetime_obj.merger_model.cost(tnode, x, multiplicity), is_log=True)
    #print(L.peak_val)
    return L, L1

# def get_old_from_new(L, node, treetime_obj, Lx=True):
#     from treetime_fft.distribution import Distribution
#     multiplicity = len(node.up.clades)
#     tnode = node.time_before_present
#     x = L.x
#     negative_merger_cost = Distribution(x, treetime_obj.merger_model.cost(tnode, x, multiplicity), is_log=True)
#     L1 = multiply(x, [L, negative_merger_cost])
#     L = Distribution(x, L.__call__(x) +treetime_obj.merger_model.cost(tnode, x, multiplicity), is_log=True)
#     #print(L.peak_val)
#     return L, L1

total_issues = 0
for key in terminals_numerical_old:
    if not terminals_numerical_old[key].time_before_present == terminals_numerical_new[key].time_before_present:
        print("time error")
    if not eq(get_new_from_old(terminals_numerical_old[key].branch_length_interpolator,terminals_numerical_old[key], tt_numerical_old )[1], terminals_numerical_new[key].branch_length_interpolator):
        if abs(get_new_from_old(terminals_numerical_old[key].branch_length_interpolator,terminals_numerical_old[key], tt_numerical_old )[1]._ymax -terminals_numerical_new[key].branch_length_interpolator._ymax) > 10e-11:
            print("_ymax too different")
            total_issues += 1
        if abs(get_new_from_old(terminals_numerical_old[key].branch_length_interpolator,terminals_numerical_old[key], tt_numerical_old )[1]._peak_val -terminals_numerical_new[key].branch_length_interpolator._peak_val) > 10e-15:
            print("_peak_val too different")
            total_issues += 1
        if abs(get_new_from_old(terminals_numerical_old[key].branch_length_interpolator,terminals_numerical_old[key], tt_numerical_old )[1].min_width -terminals_numerical_new[key].branch_length_interpolator.min_width) > 10e-11:
            print("_ymax too different")
            total_issues += 1
print("total number of issues:" +str(total_issues))

total_issues = 0
for key in names_numerical_old:
    print("number of terminals:" + str(len(names_numerical_old[key].get_terminals())))
    if not names_numerical_old[key].time_before_present == names_numerical_new[key].time_before_present:
        print("time error")
    if not eq(get_new_from_old(names_numerical_old[key].branch_length_interpolator,names_numerical_old[key], tt_numerical_old )[1], names_numerical_new[key].branch_length_interpolator):
        ymax_diff = abs(get_new_from_old(names_numerical_old[key].branch_length_interpolator,names_numerical_old[key], tt_numerical_old )[1]._ymax -names_numerical_new[key].branch_length_interpolator._ymax)
        if ymax_diff > 10e-11:
            print("_ymax too different: " + str(ymax_diff))
            total_issues += 1
        peak_val_diff = abs(get_new_from_old(names_numerical_old[key].branch_length_interpolator,names_numerical_old[key], tt_numerical_old )[1]._peak_val -names_numerical_new[key].branch_length_interpolator._peak_val)
        if peak_val_diff > 10e-15:
            print("_peak_val too different:" + str(peak_val_diff))
            total_issues += 1
        min_width_diff = abs(get_new_from_old(names_numerical_old[key].branch_length_interpolator,names_numerical_old[key], tt_numerical_old )[1].min_width -names_numerical_new[key].branch_length_interpolator.min_width)
        if min_width_diff > 10e-11:
            print("_ymax too different:" + str(min_width_diff))
            total_issues += 1
        fwhm_diff = abs(get_new_from_old(names_numerical_old[key].branch_length_interpolator,names_numerical_old[key], tt_numerical_old )[1]._fwhm -names_numerical_new[key].branch_length_interpolator._fwhm)
        if min_width_diff > 10e-15:
            print("_fwhmtoo different:" + str(fwhm_diff))
            total_issues += 1
print("total number of issues:" +str(total_issues))

# total_issues = 0
# for key in names_numerical_old:
#     print(len(names_numerical_old[key].get_terminals()))
#     if not names_numerical_old[key].time_before_present == names_numerical_new[key].time_before_present:
#         print("time error")
#     if not eq(get_old_from_new(names_numerical_new[key].branch_length_interpolator, names_numerical_new[key], tt_numerical_new )[0], names_numerical_old[key].branch_length_interpolator):
#         ymax_diff = abs(get_old_from_new(names_numerical_new[key].branch_length_interpolator,names_numerical_new[key], tt_numerical_new )[0]._ymax -names_numerical_old[key].branch_length_interpolator._ymax)
#         if ymax_diff > 10e-11:
#             print("_ymax too different: " + str(ymax_diff))
#             total_issues += 1
#         peak_val_diff = abs(get_old_from_new(names_numerical_new[key].branch_length_interpolator,names_numerical_new[key], tt_numerical_new  )[0]._peak_val -names_numerical_old[key].branch_length_interpolator._peak_val)
#         if peak_val_diff > 10e-15:
#             print("_peak_val too different:" + str(peak_val_diff))
#             total_issues += 1
#         min_width_diff = abs(get_old_from_new(names_numerical_new[key].branch_length_interpolator,names_numerical_new[key], tt_numerical_new  )[0].min_width -names_numerical_old[key].branch_length_interpolator.min_width)
#         if min_width_diff > 10e-11:
#             print("_ymax too different:" + str(min_width_diff))
#             total_issues += 1
#         fwhm_diff = abs(get_old_from_new(names_numerical_new[key].branch_length_interpolator,names_numerical_new[key], tt_numerical_new )[0]._fwhm -names_numerical_old[key].branch_length_interpolator._fwhm)
#         if min_width_diff > 10e-15:
#             print("_fwhmtoo different:" + str(fwhm_diff))
#             total_issues += 1
# print("total number of issues:" +str(total_issues))



#get the inner nodes (i.e. nodes just above the leaves)
inner_nodes = []
for key in names_numerical_old:
    if (len(names_numerical_old[key].get_nonterminals())<=1):
        inner_nodes = inner_nodes + [key]




####################This is a bunch of code for plotting the final output distributions when running both the old and the new version

tt_numerical_old = TreeTime(gtr='Jukes-Cantor', tree = tree,
                aln = aln, verbose = 1, dates = dates, root=root, seq_len =seq_len, debug=True, use_fft=False, fixed_clock_rate=fixed_clock_rate)

tt_numerical_old.run(root="best", branch_length_mode='marginal', time_marginal=True, max_iter=1, Tc=Tc)


from treetime_fft import TreeTime as TreeTimeFFT
tt_numerical_new = TreeTimeFFT(gtr='Jukes-Cantor', tree = tree,
                aln = aln, verbose = 1, dates = dates, root=root, seq_len =seq_len, debug=True, use_fft=False, fixed_clock_rate=fixed_clock_rate)

tt_numerical_new.run(root="best", branch_length_mode='marginal', time_marginal=True, max_iter=1, Tc=Tc)


def get_x_y(node, treetime_obj, Lx, new, flipped=False):
    from treetime_fft.distribution import Distribution

    if Lx:
        L = node.marginal_pos_Lx
    else:
        L = node.subtree_distribution
        x = L.x
        if new and flipped:
            multiplicity = len(node.up.clades)
            L = multiply(L.x, [L, Distribution(x, -treetime_obj.merger_model.integral_merger_rate(x), is_log=True)])
            x = L.x
    x = L.x
    if new and Lx and flipped:
        multiplicity = len(node.up.clades)
        L = multiply(L.x, [L, Distribution(x, treetime_obj.merger_model.integral_merger_rate(x), is_log=True)])
        x = L.x
    else:
        if new and Lx:
            multiplicity = len(node.up.clades)
            L = multiply(L.x, [L, Distribution(x, treetime_obj.merger_model.integral_merger_rate(x), is_log=True), Distribution(x, treetime_obj.merger_model.total_merger_rate(x)**((multiplicity-1)/(multiplicity)), is_log=False)])
            x = L.x
        else:
            if new:
                L = multiply(L.x, [L, Distribution(x, treetime_obj.merger_model.integral_merger_rate(x), is_log=True)])
                x = L.x 
    #y = np.exp(-L.y+L.peak_val)
    print(L.peak_val)
    return x, L.y

def get_new_from_old(node, treetime_obj, Lx=True):
    from treetime.distribution import Distribution
    if Lx:
        L = node.marginal_pos_Lx
        multiplicity = len(node.up.clades)
        x = L.x
        L = multiply(x, [L, Distribution(x, -treetime_obj.merger_model.integral_merger_rate(x), is_log=True), Distribution(x, treetime_obj.merger_model.total_merger_rate(x)**((-multiplicity+1)/(multiplicity)), is_log=False)])
        #y = np.exp(-L.y+L.peak_val)
    else:
        L = node.subtree_distribution
        x = L.x
        L = multiply(x, [L, Distribution(x, -treetime_obj.merger_model.integral_merger_rate(x), is_log=True)])
    print(L.peak_val)
    return x, L.y


terminals = [key for key in terminals_numerical_old]
##first check the terminal nodes, these only have the Lx (C') distribution and not the H distributions (subtree_distribution object)
pp = PdfPages('coalescent_k.pdf')
node_number = 0
for key in terminals_numerical_old:
    node_number +=1 
    print(node_number)
    if (node_number >20):
        break

    # ##calculate the x and y values for plottingcd
    x_num_1, y_num_1 = get_x_y(terminals_numerical_old[key], tt_numerical_old, True, False)
    print("number of points in old Lx:" + str(len(terminals_numerical_old[key].marginal_pos_Lx.x)))
    print("number of points in new Lx:" + str(len(terminals_numerical_new[key].marginal_pos_Lx.x)))
    x_num_3, y_num_3 = get_x_y(terminals_numerical_new[key], tt_numerical_new, True, True)
    print(x_num_1==x_num_3)
    x_num_f, y_num_f = get_x_y(terminals_numerical_new[key], tt_numerical_new, True, True, True)
    x_num_new, y_num_new = get_x_y(terminals_numerical_new[key], tt_numerical_new, True, False)
    #num_size = "numeric ultra fine, N: " + str(len(x_num_3)) + " N^2: " + str(round(len(x_num_3)**2))

    x_reverse_old, y_reverse_old = get_new_from_old(terminals_numerical_old[key], tt_numerical_old)
    x_new, y_new = get_x_y(terminals_numerical_new[key], tt_numerical_new, True, False)

    ##plot the distributions
    fig = plt.figure()
    plt.plot(x_num_1, y_num_1, 'o', markersize=1, color='red', label='numerical (old))')
    plt.plot(x_num_3, y_num_3, 'o', markersize=0.5, color='blue', label='numerical (new))')
    #plt.plot(x_num_f, y_num_f, 'o', markersize=0.5, color='black', label='numerical (flipped)')
    #plt.plot(x_num_new, y_num_new, 'o', markersize=1, color='green', label='numerical (new, not altered)')
    plt.xlabel("Time (effective support)", fontsize=8)
    plt.xlim((0, terminals_numerical_old[key].marginal_pos_Lx.peak_pos +2))
    plt.ylabel("P", fontsize=8)
    plt.title("Marginal Likelihood (Time): " + str(terminals_numerical_old[key].name), fontsize=10)
    plt.legend(prop={"size":8})
    plt.show()
    pp.savefig()
    plt.close()

    fig = plt.figure()
    plt.plot(y_num_3-y_num_1, 'o', markersize=1, color='red', label='true difference Lx old approach - new approach')
    from treetime_fft.distribution import Distribution
    integral_component = Distribution(x_num_1, tt_numerical_new.merger_model.integral_merger_rate(x_num_1), is_log=True)
    y = np.exp(-integral_component.y+integral_component.peak_val)
    multiplicity = len(terminals_numerical_new[key].up.clades)
    merger_rate_component = Distribution(x_num_1, tt_numerical_new.merger_model.total_merger_rate(x_num_1)**((multiplicity-1)/(multiplicity)), is_log=False)
    y_m = np.exp(-merger_rate_component.y+merger_rate_component.peak_val)
    plt.plot(np.max(y_num_3-y_num_1)/np.max(y)*y, 'o', markersize=1, color='green', label='integral component shape')
    plt.plot(np.max(y_num_3-y_num_1)/np.max(y_m)*y_m, 'o', markersize=1, color='blue', label='merger rate component shape')
    plt.legend(prop={"size":8})
    plt.show()
    pp.savefig()
    plt.close()

inner_nodes = []
for key in names_numerical_old:
    if (len(names_numerical_old[key].get_nonterminals())<=1):
        inner_nodes = inner_nodes + [key]

node_number = 0
for key in names_numerical_old:
    node_number +=1 
    print(node_number)
    if (node_number >20):
        break
    print(len(names_numerical_old[key].get_nonterminals()))
    if (len(names_numerical_old[key].get_nonterminals())<=1):

        ##calculate the x and y values for plottingcd
        x_num_1, y_num_1 = get_x_y(names_numerical_old[key], tt_numerical_old, False, False)
        num_size_1 = "numeric default, N: " + str(len(x_num_1)) + " N^2: " + str(round(len(x_num_1)**2))
        x_num_3, y_num_3 = get_x_y(names_numerical_new[key], tt_numerical_new, False, True)
        print("number of points in old Lx:" + str(len(names_numerical_old[key].subtree_distribution.x)))
        print("number of points in new Lx:" + str(len(names_numerical_new[key].subtree_distribution.x)))
        #x_num_f, y_num_f = get_x_y(names_numerical_new[key], tt_numerical_new, False, True, True)
        #x_num_new, y_num_new = get_x_y(names_numerical_new[key], tt_numerical_new, False, False)
        #num_size = "numeric ultra fine, N: " + str(len(x_num_3)) + " N^2: " + str(round(len(x_num_3)**2))
        x_old, y_old = get_new_from_old(names_numerical_old[key], tt_numerical_old, False)
        x_new, y_new = get_x_y(names_numerical_new[key], tt_numerical_new, False, False)
        print(y_old==y_new)

        ##plot the distributions
        fig = plt.figure()
        plt.plot(x_num_1, y_num_1, 'o', markersize=1, color='red', label='numerical (old))')
        plt.plot(x_num_3, y_num_3, 'o', markersize=0.2, color='blue', label='numerical (new))')
        #plt.plot(x_num_f, y_num_f, 'o', markersize=0.2, color='black', label='numerical (new, flipped)')
        #plt.plot(x_num_new, y_num_new, 'o', markersize=0.2, color='green', label='numerical (new, not altered)')
        plt.xlabel("Time (effective support)", fontsize=8)
        plt.xlim((0, names_numerical_old[key].marginal_pos_LH.peak_pos+10))
        plt.ylabel("P", fontsize=8)
        plt.title("Marginal Likelihood (Time): " + str(names_numerical_old[key].name), fontsize=10)
        plt.legend(prop={"size":8})
        plt.show()
        plt.savefig("checking_Lx.png")
        pp.savefig()
        plt.close()

        # ##calculate the x and y values for plottingcd
        x_num_1, y_num_1 = get_x_y(names_numerical_old[key], tt_numerical_old, True, False)
        num_size_1 = "numeric default, N: " + str(len(x_num_1)) + " N^2: " + str(round(len(x_num_1)**2))
        x_num_3, y_num_3 = get_x_y(names_numerical_new[key], tt_numerical_new, True, True)
        x_num_f, y_num_f = get_x_y(names_numerical_new[key], tt_numerical_new, True, True, True)
        x_num_new, y_num_new = get_x_y(names_numerical_new[key], tt_numerical_new, True, False)
        #num_size = "numeric ultra fine, N: " + str(len(x_num_3)) + " N^2: " + str(round(len(x_num_3)**2))

        ##plot the distributions
        fig = plt.figure()
        plt.plot(x_num_1, y_num_1, 'o', markersize=1, color='red', label='numerical (old))')
        plt.plot(x_num_3, y_num_3, 'o', markersize=1, color='blue', label='numerical (new))')
        #plt.plot(x_num_f, y_num_f, 'o', markersize=1, color='black', label='numerical (flipped)')
        plt.plot(x_num_new, y_num_new, 'o', markersize=1, color='green', label='numerical (new, not altered)')
        plt.xlabel("Time (effective support)", fontsize=8)
        plt.xlim((0, names_numerical_old[key].marginal_pos_Lx.peak_pos +10))
        plt.ylabel("P", fontsize=8)
        plt.title("Marginal Likelihood (Time): " + str(names_numerical_old[key].name), fontsize=10)
        plt.legend(prop={"size":8})
        plt.show()
        plt.savefig("checking_Lx.png")
        pp.savefig()
        plt.close()

print("done")
pp.close()

##create pre-multiplication distributions to reenact 
key=inner_nodes[0]
msgs_to_multiply_new = [names_numerical_new[key].date_constraint] if names_numerical_new[key].date_constraint is not None else []
msgs_to_multiply_new.extend([child.marginal_pos_Lx for child in names_numerical_new[key].clades
                                             if child.marginal_pos_Lx is not None])
time_points = np.unique(np.concatenate([child.marginal_pos_Lx.x for child in names_numerical_new[key].clades]))
msgs_to_multiply_new.append(tt_numerical_new.merger_model.node_contribution(names_numerical_new[key], time_points))



msgs_to_multiply_new_converted = [names_numerical_new[key].date_constraint] if names_numerical_new[key].date_constraint is not None else []
for child in names_numerical_new[key].clades:
    if child.marginal_pos_Lx is not None:
        x, y= get_x_y(child, tt_numerical_new, True, True)
        msgs_to_multiply_new_converted.extend([Distribution(x, y, is_log=True)])
msgs_to_multiply_old = [names_numerical_old[key].date_constraint] if names_numerical_old[key].date_constraint is not None else []
msgs_to_multiply_old.extend([child.marginal_pos_Lx for child in names_numerical_old[key].clades
                                             if child.marginal_pos_Lx is not None])
if(np.unique(np.concatenate([k.x for k in msgs_to_multiply_new])) == np.unique(np.concatenate([k.x for k in msgs_to_multiply_old]))):
    x_vals = np.unique(np.concatenate([k.x for k in msgs_to_multiply_old]))

m_old = multiply(x_vals, msgs_to_multiply_old )
m_new = multiply(x_vals, msgs_to_multiply_new )
from treetime.distribution import Distribution
L = multiply(x_vals, [m_old, Distribution(x_vals, -tt_numerical_new.merger_model.integral_merger_rate(x_vals), is_log=True)])