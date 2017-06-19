#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
from treetime import TreeTime, GTR
from Bio import Phylo, AlignIO

if __name__=="__main__":
    ###########################################################################
    ### parameter parsing
    ###########################################################################
    import argparse
    parser = argparse.ArgumentParser(
            description="Reconstructs ancestral sequences and infers a molecular clock tree."
                        ' The output consists of a file ending with _ancestral.fasta with ancestral sequences'
                        ' and a tree ending with _mutation.newick with mutations appended to node names'
                        ' as _A45G_.... The branches of this tree are scaled such that branch length'
                        ' correspond to times in units of the molecular clock. The molecular clock,'
                        ' along with the inferred GTR model, is written to stdout')
    parser.add_argument('--aln', required = True, type = str,  help ="fasta file with input sequences")
    parser.add_argument('--tree', required = True, type = str,  help ="newick file with tree")
    parser.add_argument('--dates', required = True, type = str,
                        help ="csv with dates for nodes with 'node_name, date' where date is float (as in 2012.15)")
    # parser.add_argument('--infer_gtr', default = True, action='store_true', help='infer substitution model')
    parser.add_argument('--gtr', required=True, type = str, help="GTR model to use. "
        " Type 'infer' to infer the model from the data. Or, specify the model type. "
        "Optionally, feed the arguments with the '--gtr_args' option")
    parser.add_argument('--gtr_params', type=str, nargs='+', help="GTR parameters for the model "
        "specified by the --gtr argument. The parameters should be feed as 'key=value' list of parameters. "
        "Example: '--gtr K80 --gtr_params kappa=0.2 pis=0.25,0.25,0.25,0.25'. See the exact definitions of "
        " the parameters in the GTR creation methods.")

    parser.add_argument('--reroot', required = False, type = str, default='best',
                        help ="reroot the tree. Valid arguments are 'best', 'midpoint', or a node name")
    parser.add_argument('--resolve_polytomies', default = True, action='store_true',
                        help='resolve polytomies using temporal information.')
    parser.add_argument('--relax',nargs='*', default = False,
                        help='use an autocorrelated molecular clock. Prior strength and coupling of parent '
                             'and offspring rates can be specified e.g. as --relax 1.0 0.5')
    parser.add_argument('--max_iter', default = 2, type=int,
                        help='maximal number of iterations the inference cycle is run')
    parser.add_argument('--verbose', default = 3, type=int,
                        help='verbosity of output 0-6')
    parser.add_argument('--Tc', default = 0.0, type=float,
                        help='coalescent time scale -- sensible values are on the order of the average '
                             'hamming distance of contemporaneous sequences')
    parser.add_argument('--plot', default = False, action='store_true',
                        help='plot the tree on a time axis')
    params = parser.parse_args()
    if params.relax==[]:
        params.relax=True
    if params.Tc<1e-5: params.Tc=False

    ###########################################################################
    ### PARSING DATES
    ###########################################################################
    with open(params.dates) as date_file:
        dates = {}
        for line in date_file:
            try:
                name, date = line.strip().split(',')[:2]
                dates[name] = float(date)
            except:
                continue


    ###########################################################################
    ### GTR SET-UP
    ###########################################################################
    model = params.gtr
    gtr_params = params.gtr_params
    if model == 'infer':
        gtr = GTR.standard('jc')
        infer_gtr = True
    else:
        kwargs = {}
        for param in gtr_params:
            keyval = param.split('=')
            if len(keyval)!=2: continue
            if keyval[0] in ['pis', 'pi', 'Pi', 'Pis']:
                keyval[1] = map(int, keyval[1].split(','))
            kwargs[keyval[0]] = keyval[1]

        gtr = GTR.standard(model, **kwargs)
        infer_gtr = False

    ###########################################################################
    ### ANCESTRAL RECONSTRUCTION AND SET-UP
    ###########################################################################
    myTree = TreeTime(dates=dates, tree=params.tree,
                       aln=params.aln, gtr=gtr, verbose=params.verbose)
    myTree.run(root=params.reroot, relaxed_clock=params.relax,
               resolve_polytomies=params.resolve_polytomies,
               Tc=params.Tc, max_iter=params.max_iter, infer_gtr=infer_gtr)

    ###########################################################################
    ### OUTPUT and saving of results
    ###########################################################################
    if params.infer_gtr:
        print('\nInferred GTR model:')
        print(myTree.gtr)

    print(myTree.date2dist)
    base_name = '.'.join(params.aln.split('/')[-1].split('.')[:-1])

    # plot
    if params.plot:
        from treetime.io import plot_vs_years
        import matplotlib.pyplot as plt
        plt.ion()
        leaf_count = myTree.tree.count_terminals()
        label_func = lambda x: x.name[:20] if (leaf_count<30 & x.is_terminal()) else ''
        branch_label_func = lambda x: (','.join([a+str(pos)+d for a,pos, d in x.mutations[:10]])
                                       +('...' if  len(x.mutations)>10 else '')) if leaf_count<30 else ''
        plot_vs_years(myTree, show_confidence=False, label_func = label_func) #, branch_labels=branch_label_func)
        plt.savefig(base_name+'_tree.pdf')
    else:
        # convert branch length to years (this is implicit in the above plot)
        myTree.branch_length_to_years()

    # decorate tree with inferred mutations
    outaln_name = base_name+'_ancestral.fasta'
    AlignIO.write(myTree.get_reconstructed_alignment(), outaln_name, 'fasta')
    for n in myTree.tree.find_clades():
        if n.up is None:
            continue
        if len(n.mutations):
            n.name+='_'+'_'.join([a+str(pos)+d for (a,pos, d) in n.mutations])

    # write tree to file. Branch length will now be scaled such that node
    # positions correspond to sampling times.
    outtree_name = base_name+'_timetree.newick'
    Phylo.write(myTree.tree, outtree_name, 'newick')

