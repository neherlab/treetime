from __future__ import print_function, division
import numpy as np
from treetime.treetime import TreeAnc as ta
from treetime.treetime import TreeTime as tt
from treetime.gtr import  GTR
from treetime import io
import datetime
import os,sys,copy
from Bio import Phylo, AlignIO
import matplotlib.pyplot as plt
plt.ion()

polytomies = False

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
    for fmt in ["%m.%d.%Y", "%m.%Y", "%Y"]:
        try:
            date = datetime.datetime.strptime(instr, fmt)
        except ValueError:
            date = None
        if date is not None:
            break
    return date

def date_from_seq_name(name):

    date = str2date_time(name.split('|')[2].strip())
    return date.year + date.timetuple().tm_yday / 365.25


if __name__=='__main__':
    root_dir = os.path.dirname(os.path.realpath(__file__))
    file_base = '../data/H3N2_NA_allyears_NA.200'
    fasta = os.path.join(root_dir, file_base+'.fasta')
    nwk = os.path.join(root_dir, file_base+'.nwk')
    mdf = os.path.join(root_dir, file_base+'.metadata.csv')

    # read tree from file
    gtr = GTR.standard()
    t = io.treetime_from_newick(gtr, nwk)
    # set alignment to the tree
    io.set_seqs_to_leaves(t, AlignIO.read(fasta, 'fasta'))
    io.read_metadata(t, mdf)
    t.reroot_to_best_root(infer_gtr=True)
    t.init_date_constraints()
    t.ml_t()
    # plotting the results
    t._score_branches()
    t.tree.ladderize()
    Phylo.draw(t.tree, label_func = lambda x:'', show_confidence=False)
    plt.title("Tree where zero-length branches are collapsed into polytomies")

    if polytomies:
        t1 = io.treetime_from_newick(gtr, nwk)
        # set alignment to the tree
        io.set_seqs_to_leaves(t1, AlignIO.read(fasta, 'fasta'))
        io.read_metadata(t1, mdf)
        t1.reroot_to_best_root(infer_gtr=True)
        t1.init_date_constraints()
        t1.ml_t()
        t1.resolve_polytomies()

        Phylo.draw(t1.tree, label_func = lambda x:'', show_confidence=False)
        plt.title("Tree with polytomies resolved")

        t.print_lh()
        print ("Prior branch len: {0}".format((t.tree.total_branch_length())))
        t1.print_lh()
        print ("Branch len after resolution: {0}".format((t1.tree.total_branch_length())))

        #traveling_wave(t1.tree, Tc=0.005)
        #t1.init_date_constraints(gtr, slope=slope)
        #t1.ml_t(gtr)
        t1.coalescent_model(optimize_Tc=True)
        t1.print_lh()
        print ("coalescent model branch len: {0}".format((t1.tree.total_branch_length())))

        gtr = GTR.standard()
        t2 = io.treetime_from_newick(gtr, nwk)
        # set alignment to the tree
        io.set_seqs_to_leaves(t2, AlignIO.read(fasta, 'fasta'))
        io.read_metadata(t2, mdf)
        t2.reroot_to_best_root(infer_gtr=True)
        t2.init_date_constraints()
        t2.ml_t()
        t2.tree.ladderize()
        t2.relaxed_clock(slack=.1, coupling=1)
        t2.ml_t()

        from matplotlib.cm import jet as cmap
        for n in t2.tree.find_clades():
            n.color = [int(x*255) for x in cmap(max(0, min(0.5*n.gamma, 1.0)))[:3]]

        Phylo.draw(t2.tree, label_func = lambda x:'', show_confidence=False, branch_labels='')
