"""
Script shows basic functionality of reading the tree and performing basic preparatory operations.

**NOTE** this script shows the  operations of the TreeAnc class, which only purpose is to prepare the tree for the ML date inferrence. The latter is done by the TreeTime class. If you want to see the functionality of this latter class, please refer to the function read_t_tree.
"""
from __future__ import print_function, division
import numpy as np
from Bio import AlignIO, Phylo
from time_tree import tree_anc as ta
import os

resources_dir = os.path.join(os.path.dirname(__file__), '../data/')

if __name__=='__main__':

    # we need the tree file and the alignment file
    tinf = os.path.join(resources_dir, 'PR.B.100.nwk')
    ainf = os.path.join(resources_dir, 'PR.B.100.fasta')

    # read alignment
    aln = AlignIO.read(ainf, 'fasta')

    # read tree from file:
    t = ta.TreeAnc.from_file(tinf, 'newick')

    # alternatively, we could construct the object directly:
    #  this is, however, limited to the valid Phylo.read  formats.
    tree_ = Phylo.read(tinf, 'newick')
    t_ = ta.TreeAnc(tree_) # Just example, not used in the following


    # first step is always to set sequences to leaves:
    err = t.set_seqs_to_leaves(aln) # err - number of failed leaves


    # As we rely on the tree provided by users, it can be built using any method/programming tool/etc., we cannot rely on the branch lengths in the input tree. they can be set in years, number of mutations, Hamming dist and what not. For us, only topology is relevant, so we want to unify all branch lengths.

    # To do this, we first make (rough) ancestral state reconstruction with Fitch algorithm (requires only topology)
    """
    **WARNING** no processing of the unknown symbols is added! If you have unknown characters, you better use ML reconstruction, which sets the unknown state to its parent state, which is always from the alphabet
    """
    t.reconstruct_anc(method='fitch')


    # When the approximate states of the internal nodes are known, we can optimize branch lengths:
    gtr = ta.GTR.standard(model='Jukes-Cantor') # model is required for ML optimization
    t.optimize_branch_len(gtr, verbose=10, store_old=True)

    # With branch lengths set in the arbitrary units, defined by the dimension of gtr.mu parameter, we can make the ML ancestral state reconstruction:
    t.reconstruct_anc(method='ml', model=gtr, verbose=10, store_lh=True)


    print ("\n\n\n\nNote that the above manipulations added the following features to the tree:")
    print ("Printing values for the root:")
    root = t.tree.root
    print ("Link to the parent node: " + str(root.up))
    print ("Distance to the root : %.3f" % root.dist2root) # depth of the node
    print ("Initial branch length (before the optimization): %.6f" %root.clades[0].old_length)
    print ("node sequence: %s " % str(root.sequence)) # sequence as numpy array
    print ("most likely state of the node (given the parent state as index) %s  " % str(root.c))
    print ("Likelihood of the most-likely node state: %s " % str(root.lh))

    # in principle, now the tree is prepared to further processing of the datetime information, which is done by the inherited TimeTree class