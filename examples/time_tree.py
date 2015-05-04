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

import sys
sys.path += ['/home/pavel/university/projects/time_tree']

import numpy as np
from Bio import AlignIO, Phylo
from time_tree import tree_anc as ta
from time_tree import time_tree as tt
import os
import datetime
resources_dir = os.path.join(os.path.dirname(__file__), '../data/')

def str2date_time(instr):
        """
        Convert input string to datetime object.

        Args:
         - instr (str): input string. Accepts one of the formats:
         {MM.DD.YYYY, MM.YYYY, MM/DD/YYYY, MM/YYYY, YYYY}.

        Returns:
         - date (datetime.datetime): parsed date object. If the parsing failed, None is returned
        """

        instr = instr.replace('/', '.')
        #import ipdb; ipdb.set_trace()
        try:
            date  = datetime.datetime.strptime(instr,  "%m.%d.%Y")
        except ValueError:
            date = None
        if date is not None:
            return date

        try:
            date  = datetime.datetime.strptime(instr,  "%m.%Y")
        except ValueError:
            date = None

        if date is not None:
            return date

        try:
            date  = datetime.datetime.strptime(instr,  "%Y")
        except ValueError:
            date = None
        return date

def flu_fasta_to_dates():
    """
    Convert fasta file with the flu data into the name,date input csv file.
    Applicable for this given format of the annotation.
    """
    ainf = os.path.join(resources_dir, 'flu.HA.fasta')  # input fasta alignment
    dinf = os.path.join(resources_dir, 'flu.HA.yrs')  # csv dates output

    outstr = []
    aln = AlignIO.read(ainf, 'fasta')
    for a in aln:
        dt = str2date_time(a.name.split('|')[2].strip())
        if dt is not None:
            outstr.append(a.name + ',' +
                    datetime.datetime.strftime(dt, "%Y.%m.%d"))

        with open (dinf, 'w') as outf:
            outf.write('\n'.join(outstr))


if __name__ == '__main__':

    # required files:
    tinf = os.path.join(resources_dir, 'flu.HA.nwk')  # newick tree
    ainf = os.path.join(resources_dir, 'flu.HA.fasta')  # fasta alignment
    dinf = os.path.join(resources_dir, 'flu.HA.yrs')  # csv dates

    # first, we need to read the tree from files.
    # There is a shortcut function, which loads data from all three files at once:
    t = tt.TreeTime.from_files(tinf, ainf, dinf, 'newick', 'fasta')

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

    # get conversion between dates and
    t.init_date_constraints(gtr)
    k_old = []
    for n in t.tree.find_clades(order='postorder'):
        k_old.append((n.dist2root, n.abs_t))
    k_old = np.array(k_old)

    # main method, performs the optimization with time constraints of the nodes
    t.ml_t(gtr)

    k_new = []
    dates = []
    for n in t.tree.find_clades(order='postorder'):
        k_new.append((n.dist2root, n.abs_t))
        dates.append((n.dist2root, n.date))
    k_new = np.array(k_new)
    dates = np.array(dates)