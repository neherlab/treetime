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
                        ' The output consists of a file ending with _ancestral.fasta with ancestral sequences'
                        ' and a tree ending with _mutation.nexus with mutations added as comments'
                        ' like _A45G_..., number in SNPs used 1-based index by default.'
                        ' The inferred GTR model is written to stdout')
    parser.add_argument('--aln', required = True, type = str,  help ="fasta file with input sequences")
    parser.add_argument('--tree', type = str,  help ="newick file with tree, "
                                                     "will attempt to build tree if none given.")

    parser.add_argument('--gtr', type = str, default='infer', help="GTR model to use. "
        " Type 'infer' to infer the model from the data. Or, specify the model type. "
        " If the specified model requires additional options, use '--gtr_args' to specify those")

    parser.add_argument('--gtr_params', type=str, nargs='+', help="GTR parameters for the model "
        "specified by the --gtr argument. The parameters should be feed as 'key=value' list of parameters. "
        "Example: '--gtr K80 --gtr_params kappa=0.2 pis=0.25,0.25,0.25,0.25'. See the exact definitions of "
        " the parameters in the GTR creation methods in treetime/nuc_models.py or treetime/aa_models.py")

    parser.add_argument('--prot', default = False, action="store_true", help ="protein alignment")
    parser.add_argument('--marginal', default = False, action="store_true", help ="marginal reconstruction of ancestral sequences")
    parser.add_argument('--zero_based', default = False, action='store_true', help='zero based SNP indexing')
    parser.add_argument('--keep_overhangs', default = False, action='store_true', help='do not fill terminal gaps')
    parser.add_argument('--verbose', default = 1, type=int, help='verbosity of output 0-6')
    params = parser.parse_args()


    ###########################################################################
    ### CHECK FOR TREE, build if not in place
    ###########################################################################
    if params.tree is None:
        from treetime.utils import tree_inference
        import os,shutil
        params.tree = os.path.basename(params.aln)+'.nwk'
        print("No tree given: inferring tree")
        tmp_dir = 'ancestral_reconstruction_tmp_files'
        tree_inference(params.aln, params.tree, tmp_dir = tmp_dir)
        if os.path.isdir(tmp_dir):
            shutil.rmtree(tmp_dir)

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
    treeanc = TreeAnc(params.tree, aln=params.aln, gtr=gtr, verbose=4, fill_overhangs=not params.keep_overhangs)
    treeanc.infer_ancestral_sequences('ml', infer_gtr=infer_gtr,
                                       marginal=params.marginal)

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
