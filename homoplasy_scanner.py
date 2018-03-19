#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
from treetime import TreeAnc, GTR
from Bio import Phylo, AlignIO
from Bio import __version__ as bioversion

if __name__=="__main__":
    ###########################################################################
    ### parameter parsing
    ###########################################################################
    import argparse
    parser = argparse.ArgumentParser(
            description='Reconstructs ancestral sequences and maps mutations to the tree.'
                        ' The tree is then scanned for homoplasies. An excess number of homoplasies'
                        ' might suggest contamination, recombination, culture adaptation or similar. ')
    parser.add_argument('--aln', required = True, type = str,  help ="fasta file with input sequences")
    parser.add_argument('--tree', required = True, type = str,  help ="newick file with tree")

    parser.add_argument('--gtr', required=False, type = str, default='infer', help="GTR model to use. "
        " Type 'infer' to infer the model from the data. Or, specify the model type. "
        " If the specified model requires additional options, use '--gtr_args' to specify those")

    parser.add_argument('--gtr_params', type=str, nargs='+', help="GTR parameters for the model "
        "specified by the --gtr argument. The parameters should be feed as 'key=value' list of parameters. "
        "Example: '--gtr K80 --gtr_params kappa=0.2 pis=0.25,0.25,0.25,0.25'. See the exact definitions of "
        " the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py")

    parser.add_argument('--prot', default = False, action="store_true", help ="protein alignment")
    parser.add_argument('--marginal', default = False, action="store_true", help ="protein alignment")
    parser.add_argument('--zero_based', default = False, action='store_true', help='zero based SNP indexing')
    parser.add_argument('--verbose', default = 1, type=int, help='verbosity of output 0-6')
    params = parser.parse_args()



    ###########################################################################
    ### GTR SET-UP
    ###########################################################################
    model = params.gtr
    gtr_params = params.gtr_params
    if model == 'infer':
        gtr = GTR.standard('jc')
        infer_gtr = True
    else:
        try:
            kwargs = {}
            if gtr_params is not None:
                for param in gtr_params:
                    keyval = param.split('=')
                    if len(keyval)!=2: continue
                    if keyval[0] in ['pis', 'pi', 'Pi', 'Pis']:
                        keyval[1] = map(float, keyval[1].split(','))
                    elif keyval[0] not in ['alphabet']:
                        keyval[1] = float(keyval[1])
                    kwargs[keyval[0]] = keyval[1]
            else:
                print ("GTR params are not specified. Creating GTR model with default parameters")


            gtr = GTR.standard(model, **kwargs)
            infer_gtr = False
        except:
            print ("Could not create GTR model from input arguments. Using default (Jukes-Cantor 1969)")
            gtr = GTR.standard('jc')
            infer_gtr = False


    ###########################################################################
    ### ANCESTRAL RECONSTRUCTION
    ###########################################################################
    treeanc = TreeAnc(params.tree, aln=params.aln, gtr=gtr, verbose=1,
                      fill_overhangs=True)
    L = treeanc.aln.get_alignment_length()
    N_seq = len(treeanc.aln)
    N_tree = treeanc.tree.count_terminals()

    print("read alignment from file %s with %d sequences of length %d"%(params.aln,N_seq,L))
    print("read tree from file %s with %d leaves"%(params.tree,N_tree))

    treeanc.infer_ancestral_sequences('ml', infer_gtr=infer_gtr,
                                       marginal=params.marginal)

    ###########################################################################
    ### analysis of reconstruction
    ###########################################################################
    from collections import defaultdict
    from scipy.stats import poisson
    offset = 0 if params.zero_based else 1

    mutations = defaultdict(list)
    positions = defaultdict(list)
    terminal_mutations = defaultdict(list)

    for n in treeanc.tree.find_clades():
        if n.up is None:
            continue

        if len(n.mutations):
            for (a,pos, d) in n.mutations:
                if '-' not in [a,d]:
                    mutations[(a,pos+offset,d)].append(n)
                    positions[pos+offset].append(n)
            if n.is_terminal():
                for (a,pos, d) in n.mutations:
                    if '-' not in [a,d]:
                        terminal_mutations[(a,pos+offset,d)].append(n)

    total_branch_length = treeanc.tree.total_branch_length()
    corrected_branch_length = np.sum([np.exp(-x.branch_length)*np.sinh(x.branch_length)
                                      for x in treeanc.tree.find_clades()])
    corrected_terminal_branch_length = np.sum([np.exp(-x.branch_length)*np.sinh(x.branch_length)
                                      for x in treeanc.tree.get_terminals()])
    expected_mutations = L*corrected_branch_length
    expected_terminal_mutations = L*corrected_terminal_branch_length
    multiplicities = np.bincount([len(x) for x in mutations.values()])
    total_mutations = np.sum([len(x) for x in mutations.values()])

    multiplicities_terminal = np.bincount([len(x) for x in terminal_mutations.values()])
    terminal_mutation_count = np.sum([len(x) for x in terminal_mutations.values()])

    print("\nThe TOTAL tree length is %1.3e, expecting %1.1f mutations vs an observed %d"
          %(total_branch_length,expected_mutations,total_mutations))
    print("Of these %d mutations,"%total_mutations
            +"".join(['\n\t - %d occur %d times'%(n,mi)
                      for mi,n in enumerate(multiplicities) if n]))
    print("\nThe TERMINAL branch length is %1.3e, expecting %1.1f mutations vs an observed %d"
          %(corrected_terminal_branch_length,expected_terminal_mutations,terminal_mutation_count))
    print("Of these %d mutations,"%terminal_mutation_count
            +"".join(['\n\t - %d occur %d times'%(n,mi)
                      for mi,n in enumerate(multiplicities_terminal) if n]))


    multiplicities_positions = np.bincount([len(x) for x in positions.values()])
    multiplicities_positions[0] = L - np.sum(multiplicities_positions)
    print("\nOf the %d positions in the genome,"%L
            +"".join(['\n\t - %d were hit %d times (expected %1.2f)'%(n,mi,L*poisson.pmf(mi,corrected_branch_length))
                      for mi,n in enumerate(multiplicities_positions) if n]))


    p = poisson.pmf(np.arange(10*multiplicities_positions.max()),1.0*total_mutations/L)
    print("\nlog-likelihood difference to Poisson distribution with same mean: %1.3e"%(
            - L*np.sum(p*np.log(p+1e-100))
            + np.sum(multiplicities_positions*np.log(p[:len(multiplicities_positions)]+1e-100))))


    print("\n\nThe ten most homoplasic mutations are:\n\tmut\tmultiplicity")
    mutations_sorted = sorted(mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
    for mut, val in mutations_sorted[:10]:
        if len(val)>1:
            print("\t%s%d%s\t%d"%(mut[0], mut[1], mut[2], len(val)))
        else:
            break

    print("\n\nThe ten most homoplasic mutation on terminal branches are:\n\tmut\tmultiplicity")
    terminal_mutations_sorted = sorted(terminal_mutations.items(), key=lambda x:len(x[1])-0.1*x[0][1]/L, reverse=True)
    for mut, val in terminal_mutations_sorted[:10]:
        if len(val)>1:
            print("\t%s%d%s\t%d"%(mut[0], mut[1], mut[2], len(val)))
        else:
            break

