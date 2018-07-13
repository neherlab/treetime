#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
from treetime import TreeAnc, GTR
from utils import assure_tree, create_gtr
from treetime import config as ttconf
from Bio import Phylo, AlignIO
from Bio import __version__ as bioversion
import sys

def ancestral(params)
    ###########################################################################
    ### CHECK FOR TREE, build if not in place
    ###########################################################################
    if assure_tree(params.tree, tmp_dir='homoplasy_tmp'):
        return 1

    ###########################################################################
    ### GTR SET-UP
    ###########################################################################
    gtr = create_gtr(params)

    ###########################################################################
    ### ANCESTRAL RECONSTRUCTION
    ###########################################################################
    treeanc = TreeAnc(params.tree, aln=params.aln, gtr=gtr, verbose=4,
                      fill_overhangs=not params.keep_overhangs)
    ndiff =treeanc.infer_ancestral_sequences('ml', infer_gtr=params.gtr=='infer',
                                             marginal=params.marginal)
    if ndiff==ttconf.ERROR: # if reconstruction failed, exit
        sys.exit(1)

    ###########################################################################
    ### OUTPUT and saving of results
    ###########################################################################

    model = 'aa' if params.prot else 'Jukes-Cantor'
    if infer_gtr:
        print('\nInferred GTR model:')
        print(treeanc.gtr)

    outaln_name = '.'.join(params.aln.split('/')[-1].split('.')[:-1])+'_ancestral.fasta'
    AlignIO.write(treeanc.get_reconstructed_alignment(), outaln_name, 'fasta')
    print("--- alignment including ancestral nodes saved as  \n\t %s\n"%outaln_name)

    # decorate tree with inferred mutations
    terminal_count = 0
    offset = 0 if params.zero_based else 1
    for n in treeanc.tree.find_clades():
        if n.up is None:
            continue
        n.confidence=None
        # due to a bug in older versions of biopython that truncated filenames in nexus export
        # we truncate them by hand and make them unique.
        if n.is_terminal() and len(n.name)>40 and bioversion<"1.69":
            n.name = n.name[:35]+'_%03d'%terminal_count
            terminal_count+=1
        if len(n.mutations):
            n.comment= '&mutations="' + '_'.join([a+str(pos + offset)+d for (a,pos, d) in n.mutations])+'"'

    # write tree to file
    outtree_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])+'_mutation.nexus'
    Phylo.write(treeanc.tree, outtree_name, 'nexus')
    print("--- tree saved in nexus format as  \n\t %s\n"%outtree_name)

    return 0
