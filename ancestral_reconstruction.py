#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
from treetime import TreeAnc
from Bio import Phylo, AlignIO

if __name__=="__main__":
    ###########################################################################
    ### parameter parsing
    ###########################################################################
    import argparse
    parser = argparse.ArgumentParser(
            description='Reconstructs ancestral sequences and maps mutations to the tree.'
                        ' The output consists of a file ending with _ancestral.fasta with ancestral sequences'
                        ' and a tree ending with _mutation.nexus with mutations added as comments'
                        ' like _A45G_... The inferred GTR model is written to stdout')
    parser.add_argument('--aln', required = True, type = str,  help ="fasta file with input sequences")
    parser.add_argument('--tree', required = True, type = str,  help ="newick file with tree")
    parser.add_argument('--prot', default = False, action="store_true", help ="protein alignment")
    parser.add_argument('--marginal', default = False, action='store_true', help='marginal instead of joint ML reconstruction')
    parser.add_argument('--infer_gtr', default = False, action='store_true', help='infer substitution model')
    parser.add_argument('--keep_overhangs', default = False, action='store_true', help='do not fill terminal gaps')
    params = parser.parse_args()

    ###########################################################################
    ### ANCESTRAL RECONSTRUCTION
    ###########################################################################
    model = 'aa' if params.prot else 'Jukes-Cantor'
    treeanc = TreeAnc(params.tree, aln=params.aln, gtr=model, verbose=4, fill_overhangs=not params.keep_overhangs)
    treeanc.infer_ancestral_sequences('ml', infer_gtr=params.infer_gtr,
                                       marginal=params.marginal)

    ###########################################################################
    ### OUTPUT and saving of results
    ###########################################################################
    if params.infer_gtr:
        print('\nInferred GTR model:')
        print(treeanc.gtr)

    outaln_name = '.'.join(params.aln.split('/')[-1].split('.')[:-1])+'_ancestral.fasta'
    AlignIO.write(treeanc.get_reconstructed_alignment(), outaln_name, 'fasta')

    # decorate tree with inferred mutations
    terminal_count = 0
    for n in treeanc.tree.find_clades():
        if n.up is None:
            continue
        n.confidence=None
        if n.is_terminal() and len(n.name)>40:
            n.name = n.name[:35]+'_%03d'%terminal_count
            terminal_count+=1
        if len(n.mutations):
            n.comment= '&mutations="' + '_'.join([a+str(pos)+d for (a,pos, d) in n.mutations])+'"'

    # write tree to file
    outtree_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])+'_mutation.nexus'
    Phylo.write(treeanc.tree, outtree_name, 'nexus')
