from __future__ import print_function, division
import numpy as np
from Bio import Phylo
from treetime import TreeTime
from treetime.utils import parse_dates
from treetime.node_interpolator import Distribution, NodeInterpolator


if __name__ == '__main__':

    base_name = 'test/treetime_examples/data/h3n2_na/h3n2_na_20'

    dates = parse_dates(base_name+'.metadata.csv')
    tt = TreeTime(gtr='Jukes-Cantor', tree = base_name+'.nwk',
                  aln = base_name+'.fasta', verbose = 3, dates = dates, debug=True)

    # rerooting can be done along with the tree time inference
    tt.run(root="best", branch_length_mode='input', max_iter=2, time_marginal=True)

    # initialize date constraints and branch length interpolators
    # this called in each iteration. 44ms
    tt.init_date_constraints()

    ###########################################################
    # joint inference of node times. done in every generation. 0.7s
    tt._ml_t_joint()
    # individual steps in joint inference - post-order
    msgs_to_multiply = [child.joint_pos_Lx for child in tt.tree.root.clades
                                             if child.joint_pos_Lx is not None]

    # 330us
    subtree_distribution = Distribution.multiply(msgs_to_multiply)
    # 30ms (there are 19 nodes here, so about 20 internal branches -> 1s)
    res, res_t = NodeInterpolator.convolve(subtree_distribution, tt.tree.root.clades[1].branch_length_interpolator, max_or_integral='max', inverse_time=True, n_grid_points = tt.node_grid_points, n_integral=tt.n_integral, rel_tol=tt.rel_tol_refine)


    ###########################################################
    # marginal inference. done only for confidence estimation: 2.7s
    tt._ml_t_marginal()

    # individual steps in marginal inference - post-order
    msgs_to_multiply = [child.marginal_pos_Lx for child in tt.tree.root.clades
                                             if child.marginal_pos_Lx is not None]

    # 330us
    subtree_distribution = Distribution.multiply(msgs_to_multiply)
    # 60ms (there are 19 nodes here, so about 20 internal branches -> 1s)
    res, res_t = NodeInterpolator.convolve(subtree_distribution, tt.tree.root.clades[1].branch_length_interpolator, max_or_integral='integral', inverse_time=True, n_grid_points = tt.node_grid_points, n_integral=tt.n_integral, rel_tol=tt.rel_tol_refine)
    # 80ms (there are 19 nodes here, so about 20 internal branches -> 1s)
    res, res_t = NodeInterpolator.convolve(subtree_distribution, tt.tree.root.clades[1].branch_length_interpolator, max_or_integral='integral', inverse_time=False, n_grid_points = tt.node_grid_points, n_integral=tt.n_integral, rel_tol=tt.rel_tol_refine)

    # This points towards the convolution being the biggest computational expense. 



