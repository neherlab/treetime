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

if __name__=='__main__':

    gtr = ta.GTR.standard()

    root_dir = os.path.dirname(os.path.realpath(__file__))
    fasta = os.path.join(root_dir, '../data/flu_trivial.fasta')
    nwk = os.path.join(root_dir, '../data/flu_trivial.nwk')

    # read tree from file
    t = io.treetime_from_newick(nwk)
    io.set_seqs_to_leaves(t, AlignIO.read(fasta, 'fasta'))
    io.set_node_dates_from_names(t, date_from_seq_name)

    t.reroot_to_oldest()
    t.optimize_seq_and_branch_len(gtr,reuse_branch_len=True,prune_short=True)
    
    slope = None #  -1.1505574145108622e-05
    t.init_date_constraints(gtr, slope=slope)
    
    t.ml_t(gtr)
    # plotting the results

    t._score_branches()
    Phylo.draw(t.tree, label_func=lambda x: round(x.numdate,2), show_confidence=False, branch_labels=None)
    plt.savefig("polytomies_raw.svg")

    t.resolve_polytomies(gtr, False)
    t._score_branches()
    Phylo.draw(t.tree, label_func=lambda x: round(x.numdate,2), show_confidence=False, branch_labels=None)
    plt.savefig("polytomies_stretched_only.svg")

    t.resolve_polytomies(gtr, True)
    t._score_branches()
    Phylo.draw(t.tree, label_func=lambda x: round(x.numdate,2), show_confidence=False, branch_labels=None)
    plt.savefig("polytomies_stretched_comp.svg")




