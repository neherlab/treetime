from __future__ import print_function, division

import numpy as np
from tree_time.tree_time import tree_anc as ta
from tree_time.tree_time import tree_time as tt
import datetime
import os
from Bio import Phylo, AlignIO
import matplotlib.pyplot as plt
plt.ion()


def str2date_time(instr):
    """
    Convert input string to datetime object.

    Args:
     - instr (str): input string. Accepts one of the formats:
     {MM.DD.YYYY, MM.YYYY, MM/DD/YYYY, MM/YYYY, YYYY}.

    Returns:
     - date (datetime.datetime): parsed date object. If the parsing failed,
     None is returned
    """

    instr = instr.replace('/', '.')
    # import ipdb; ipdb.set_trace()
    try:
        date = datetime.datetime.strptime(instr, "%m.%d.%Y")
    except ValueError:
        date = None
    if date is not None:
        return date

    try:
        date = datetime.datetime.strptime(instr, "%m.%Y")
    except ValueError:
        date = None

    if date is not None:
        return date

    try:
        date = datetime.datetime.strptime(instr, "%Y")
    except ValueError:
        date = None
    return date

def date_from_seq_name(name):
    
    date = str2date_time(name.split('|')[2].strip())

    return date

def binstr(x):
    if x > 100:
        return '>10 or <0'
    else:
        return '<%.2f' % x

def cost_fun(n):
    sign = np.sign(n.branch_length - t.opt_branch_length(n))
    return sign * (n.branch_neg_log_prob(t.opt_branch_length(n)) - n.branch_neg_log_prob(n.branch_length))


if __name__=='__main__':

    
    gtr = ta.GTR.standard()

    root_dir = os.path.dirname(os.path.realpath(__file__)) 
    fasta = os.path.join(root_dir, 'H3N2_NA_allyears_NA.200.fasta')
    nwk = os.path.join(root_dir, 'H3N2_NA_allyears_NA.200.nwk')

    # read tree from file
    t = tt.TreeTime.from_newick(nwk)
    
    # set alignment to the tree
    aln = AlignIO.read(fasta, 'fasta')
    t.load_aln(aln)
    # set dates from the node names 
    t.set_node_dates_from_names(date_from_seq_name)

    
    t.reroot_to_oldest()

    t.optimize_seq_and_branch_len(gtr) 
    #import ipdb; ipdb.set_trace()

    rds = []
    for n in t.tree.find_clades():
        if hasattr(n, 'raw_date') and n.raw_date is not None: 
            rds.append((n.dist2root, n.raw_date ))

    
    t.init_date_constraints(gtr, slope=None)
    t.ml_t(gtr)

    # plotting the results
    
    t._score_branches()
    #scores = []
    ##for n in t.tree.find_clades():
    ##    scores.append(n.score)
    ##
    #k = map(binstr, bins)
    
    #plt.figure(1)
    #plt.hist(scores)
    #plt.xticks(np.linspace(0,0.9,10) + 0.05, list(k))
    #plt.ylabel('Branch count')
    #plt.xlabel('Branch length deviation from optimum,\nrelated to the optimal branch length')
    #plt.title ('The deviation of the branch lengths from their optimal values\nInfluenza H3N2')

    #
    t.tree.ladderize()
    Phylo.draw(t.tree, label_func = lambda x:'', show_confidence=False, branch_labels='')
    t.print_lh()
    
    t.resolve_polytomies(gtr, t.ladderize_node_polytomies, (gtr,))
    t._score_branches()
    Phylo.draw(t.tree, label_func = lambda x:'', show_confidence=False, branch_labels='')
    t.print_lh()
    