from __future__ import print_function, division
import numpy as np
from tree_time.tree_time import tree_anc as ta
from tree_time.tree_time import tree_time as tt
from tree_time.tree_time.gtr import GTR
from tree_time.tree_time import io
from tree_time.tree_time.merger_models import coalescent
import datetime
import os,sys,copy
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

if __name__=='__main__':
    gtr = GTR.standard()
    root_dir = os.path.dirname(os.path.realpath(__file__))
    fasta = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.200.fasta')
    nwk = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.200.nwk')
    #fasta = os.path.join(root_dir, 'flu_trivial.fasta')
    #nwk = os.path.join(root_dir, 'flu_trivial.nwk')
    slope = -1.1505574145108622e-05
    # read tree from file
    t = io.treetime_from_newick(nwk)
    # set alignment to the tree
    io.set_seqs_to_leaves(t, AlignIO.read(fasta, 'fasta'))
    # set dates from the node names
    io.set_node_dates_from_names(t, date_from_seq_name)
    t.reroot_to_oldest()
    t.optimize_seq_and_branch_len(gtr)
    t.init_date_constraints(gtr, slope=slope)
    t.ml_t(gtr)
    # plotting the results
    t._score_branches()
    t.tree.ladderize()
    coalescent(t.tree, Tc=0.1)
    #Phylo.draw(t.tree, label_func = lambda x:'', show_confidence=False, branch_labels='')
    t1 = copy.deepcopy(t)
    t1.resolve_polytomies(gtr)
    t1.tree.ladderize()
    t.print_lh()
    print ("Prior branch len: {0}".format((t.tree.total_branch_length())))
    t1.print_lh()
    print ("Posterior branch len: {0}".format((t1.tree.total_branch_length())))
