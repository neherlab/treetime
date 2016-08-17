from __future__ import print_function, division
import numpy as np
from treetime.clock_tree import ClockTree
from Bio import Phylo, AlignIO

if __name__=="__main__":
    ###########################################################################
    ### parameter parsing
    ###########################################################################
    import argparse
    parser = argparse.ArgumentParser(
            description="Reconstruct ancestral sequences, set dates to tree, and infer a time scaled tree."
                        "The tree needs to be properly rooted -- other than branch length the tree won't be modified"
                        " The ancestral sequences will be written to a file ending on _ancestral.fasta"
                        " A tree in newick format with mutations as _A45G_... appended"
                        " appended to node names will be written to a file ending on _mutation.newick")
    parser.add_argument('--aln', required = True, type = str,  help ="fasta file with input sequences")
    parser.add_argument('--tree', required = True, type = str,  help ="newick file with tree")
    parser.add_argument('--dates', required = True, type = str,  help ="csv with dates (float as in 2012.15) for nodes")
    parser.add_argument('--infer_gtr', default = False, action='store_true', help='infer substitution model')
    params = parser.parse_args()

    ###########################################################################
    ### PARSING DATES
    ###########################################################################
    with open(params.dates) as date_file:
        dates = {}
        for line in date_file:
            try:
                name, date = line.strip().split(',')
                dates[name] = float(date)
            except:
                continue

    ###########################################################################
    ### ANCESTRAL RECONSTRUCTION AND SET-UP
    ###########################################################################
    myTree = ClockTree(dates, params.tree, aln=params.aln, gtr='Jukes-Cantor', verbose=4)
    myTree.init_date_constraints(infer_gtr=params.infer_gtr)

    ###########################################################################
    ### OUTPUT and saving of results
    ###########################################################################
    if params.infer_gtr:
        print('\nInferred GTR model:')
        print(myTree.gtr)

    print(myTree.date2dist)

    outaln_name = '.'.join(params.aln.split('/')[-1].split('.')[:-1])+'_ancestral.fasta'
    AlignIO.write(myTree.get_reconstructed_alignment(), outaln_name, 'fasta')

    # decorate tree with inferred mutations
    for n in myTree.tree.find_clades():
        if n.up is None:
            continue
        if len(n.mutations):
            n.name+='_'+'_'.join([a+str(pos)+d for (a,pos, d) in n.mutations])

    # write tree to file
    outtree_name = '.'.join(params.tree.split('/')[-1].split('.')[:-1])+'_mutation.newick'
    Phylo.write(myTree.tree, outtree_name, 'newick')
