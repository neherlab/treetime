from __future__ import print_function, division
import numpy as np
from treetime.treetime import treeanc as ta
from treetime.treetime import treetime as tt
from treetime.treetime.gtr import GTR
from treetime.treetime import io
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

def optimize_tree(nwk, fasta, gtr, slope=None):
    """
    Read data from text files and do the tree processing pipeline.
    Returns:
     - tree (TreeTime): optimized tree.
    """
    # initialize the tree first
    t = io.treetime_from_newick(nwk)
    io.set_seqs_to_leaves(t, AlignIO.read(fasta, 'fasta'))
    io.set_node_dates_from_names(t, date_from_seq_name)

    # prepare branch lengths and sequences of the internal nodes
    t.reroot_to_oldest()
    t.optimize_seq_and_branch_len(gtr)
    # prepare the date-time information about the nodes
    t.init_date_constraints(gtr, slope=slope)
    # do the optimization
    t.ml_t(gtr)
    t.resolve_polytomies(gtr)
    t.print_lh()
    return t

def lh_for_slope(t, slope):
    """
    Compute the tree - likelihood for a given slope
    """
    t1 = copy.deepcopy(t) # should be redundant, just to test
    t1.init_date_constraints(gtr, slope=slope)
    t1.ml_t(gtr)
    return t1.total_LH()

if __name__=='__main__':
    #fasta = os.path.join(root_dir, 'flu_trivial.fasta')
    #nwk = os.path.join(root_dir, 'flu_trivial.nwk')
    gtr = GTR.standard()
    root_dir = os.path.dirname(os.path.realpath(__file__))
    fasta = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.20.fasta')
    nwk = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.20.nwk')
    slope = -1.1505574145108622e-05
    # read tree from file
    t = optimize_tree(nwk, fasta, gtr, slope)
    slopes = slope + slope*np.linspace(-0.5, 0.5, 50) # 50% slope variation
    res = []
    i = 0
    for sl in slopes:
        i += 1
        print ("Iteration:" + str(i) + " Tree likelihood calculation for slope: " + str(sl))
        res.append(lh_for_slope(t, sl))
    plt.plot(slopes, res)
