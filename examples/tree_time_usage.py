"""
Script shows basic functionality of reading the tree and doing the ML
optimization given the datetime constraints of the nodes.
To perform the optimization, a special class TimeTree was created. This
example script shows howe to use this class.

**NOTE** this script shows only the examples for  TimeTree class, which is
supposed to do **only** the ML optimization with constraints (+ necessary
preparation). For basic functionality, reading-writing, inferring
ancestral states, see tree_anc.py example and refer to TreeAnc class (base
class of TimeTree) and corresponding example.
"""

from __future__ import print_function, division

#import sys
#sys.path += ['/home/pavel/university/projects']

#from .read_tree import W

#sys.exit(0)
import numpy as np
from Bio import AlignIO, Phylo
from tree_time.tree_time import tree_anc as ta
from tree_time.tree_time import tree_time as tt
import os
import datetime
resources_dir = os.path.join(os.path.dirname(__file__), '../data/')

def n_name(node):
    s_l = ''
    s_r = ''
    s_b = ''
    s_opt = ''
    if hasattr(node, 'branch_neg_log_prob') and node.branch_neg_log_prob is not None:
        log = node.neg_log_prob
        x = node.abs_t
        s_l = '%.0f' % ((log(x) - log(x - 1e-5)) / 1e-5)
        s_r = '%.0f' % ((log(x + 1e-5) - log(x)) / 1e-5)
        
        x1 = node.branch_length
        s_opt = '%.0f' % ((log(x1+1e-5) - log(x1-1e-5)) / 2e-5)
        log = node.branch_neg_log_prob
        s_b = '%.4f' % (log.x[log(log.x).argmin()] - node.branch_length)
    if node.up is None:
        s2 = ''
    else:
        s2 = str(node.profile.shape[0] - (node.profile * node.up.profile).sum())
    return s_opt +'//' + s_l + '//' +  s_r + '//' + s_b + '//' + s2


if __name__ == '__main__':

    # required files:
    tinf = os.path.join(resources_dir, 'flu.HA.nwk')  # newick tree
    ainf = os.path.join(resources_dir, 'flu.HA.fasta')  # fasta alignment
    dinf = os.path.join(resources_dir, 'flu.HA.yrs')  # csv dates

    # first, we need to read the tree from files.
    # There is a shortcut function, which loads data from all three files at
    # once:
    t = tt.TreeTime.from_files(tinf, ainf, dates_file=dinf)

    # we can of course do everything manually (just as an example)
    # _t = tt.TreeTime.from_file(tinf, 'newick')
    # note we call function in parent class (TreeAnc) from chld class (TreeTime)
    # _aln = AlignIO.read(ainf, 'fasta')
    # _t.set_seqs_to_leaves(_aln)  # load sequences from alignment
    # _t.load_dates(dinf)  # load dates from csv file
    # now, we need preparation. for details, see tree_anc.py example.
    # normally, we start with Fitch reconstruction,
    # but since it does not support unknown characters, we do ML:
    gtr = ta.GTR.standard(model='Jukes-Cantor')  # model is required for ML
    t.reconstruct_anc('ml', model=gtr)
    #  FIXME sometimes fails

    t.optimize_branch_len(gtr, verbose=10, store_old=True)
    # and reconsruct  once more:
    # t.reconstruct_anc('ml',model=gtr)

    # now we are ready to make the tree optimization with time constraints

    t.prune_short_branches()
    # get conversion between dates and
    t.init_date_constraints(gtr)
    k_old = []
    for n in t.tree.find_clades(order='postorder'):
        if hasattr(n, 'raw_date') and n.raw_date is not None:
            k_old.append((n.raw_date, n.dist2root))
    k_old = np.array(k_old)

    # main method, performs the optimization with time constraints of the nodes
    t.ml_t(gtr)

    k_new = []
    dates = []
    for n in t.tree.find_clades(order='postorder'):
        k_new.append((n.date, n.abs_t))
        dates.append((n.dist2root, n.date))
    k_new = np.array(k_new)
    dates = np.array(dates)
