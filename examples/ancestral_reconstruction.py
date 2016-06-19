from __future__ import print_function, division
import numpy as np
from treetime import io
from treetime.gtr import  GTR
from Bio import Phylo, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(
            description='Reconstruct ancestral sequences and map mutations to the tree.'
                        ' The ancestral sequences will be written to a file "aln_base"_ancestral.fasta'
                        ' A tree in newick format with mutations as _A45G_... appended'
                        ' appended to node names will be written to a file "treebase"_mutation.newick')
    parser.add_argument('--aln', required = True, type = str,  help ="fasta file with input sequences")
    parser.add_argument('--tree', required = True, type = str,  help ="newick file with tree")
    parser.add_argument('--marginal', default = False, action='store_true', help='marginal instead of joint ML reconstruction')
    parser.add_argument('--infer_gtr', default = False, action='store_true', help='infer substitution model')

    params = parser.parse_args()

    gtr = GTR.standard()
    treeanc = io.treetime_from_newick(gtr, params.tree)
    io.set_seqs_to_leaves(treeanc, AlignIO.read(params.aln, 'fasta'))

    treeanc.reconstruct_anc('ml', infer_gtr=params.infer_gtr, marginal=params.marginal)
    if params.infer_gtr:
        print('Inferred GTR model:')
        print(treeanc.gtr)


    new_aln = MultipleSeqAlignment([SeqRecord(id=n.name, seq=Seq("".join(n.sequence)),
                                              description="")
                                    for n in treeanc.tree.find_clades()])
    outaln_name = '.'.join(params.aln.split('/')[-1].split('.')[:-1])+'_ancestral.fasta'
    AlignIO.write(new_aln, outaln_name, 'fasta')

    for n in treeanc.tree.find_clades():
        if n.up is None:
            continue
        if len(n.mutations):
            n.name+='_'+'_'.join([a+str(pos)+d for (a,pos, d) in n.mutations])

    outtree_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])+'_mutation.newick'
    Phylo.write(treeanc.tree, outtree_name, 'newick')
