#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
import pandas as pd
from treetime import TreeAnc, GTR
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import Phylo, AlignIO
from Bio import __version__ as bioversion
import os

if __name__=="__main__":
    ###########################################################################
    ### parameter parsing
    ###########################################################################
    import argparse
    parser = argparse.ArgumentParser(
            description='Reconstructs discrete ancestral states, for example '
                        'geographic location, host, or similar.')
    parser.add_argument('--tree', required = True, type=str, help ="newick file with tree")
    parser.add_argument('--attribute', type=str, help ="attribute to reconstruct, e.g. country")
    parser.add_argument('--states', required = True, type=str, help ="csv or tsv file with discrete characters."
                                    "\n#name,country,continent\ntaxon1,micronesia,oceania\n...")
    parser.add_argument('--weights', type=str, help="csv or tsv file with probabilities of that a randomly sampled "
                        "sequence has a particular state. E.g. population of different continents or countries. E.g.:"
                        "\n#country,weight\nmicronesia,0.1\n...")
    # parser.add_argument('--migration', type=str, help="csv or tsv file with symmetric migration/transition rates "
    #                     "between states. For example passenger flow.")
    # parser.add_argument('--infer_gtr', action="store_true", help="infer GTR model from tree. "
    #                                 "Ignored when prop or migration is specified.")
    parser.add_argument('--confidence', action="store_true", help="output confidence of mugration inference")

    parser.add_argument('--verbose', default = 1, type=int, help='verbosity of output 0-6')
    params = parser.parse_args()
    missing = "?"

    ###########################################################################
    ### Parse states
    ###########################################################################
    if os.path.isfile(params.states):
        states = pd.read_csv(params.states, sep='\t' if params.states[-3:]=='tsv' else ',',
                             skipinitialspace=True)
    else:
        print("file with states does not exist")
        exit(1)

    taxon_name = 'name' if 'name' in states.columns else states.columns[0]
    if params.attribute and params.attribute in states.columns:
        attr = params.attribute
    else:
        attr = states.columns[1]

    leaf_to_attr = {x[taxon_name]:x[attr] for xi, x in states.iterrows()
                    if x[attr]!=missing}
    unique_states = sorted(set(leaf_to_attr.values()))
    nc = len(unique_states)
    if nc>180:
        print("mugration: can't have more than 180 states!")
        exit(1)
    elif nc<2:
        print("mugration: only one or zero states found -- this doesn't make any sense")
        exit(1)

    ###########################################################################
    ### make a single character alphabet that maps to discrete states
    ###########################################################################
    alphabet = [chr(65+i) for i,state in enumerate(unique_states)]
    missing_char = chr(65+nc)
    letter_to_state = {a:unique_states[i] for i,a in enumerate(alphabet)}
    letter_to_state[missing_char]=missing
    reverse_alphabet = {v:k for k,v in letter_to_state.items()}

    ###########################################################################
    ### construct gtr model
    ###########################################################################
    if params.weights:
        params.infer_gtr = True
        tmp_weights = pd.read_csv(params.weights, sep='\t' if params.states[-3:]=='tsv' else ',',
                             skipinitialspace=True)
        weights = {row[0]:row[1] for ri,row in tmp_weights.iterrows()}
        mean_weight = np.mean(weights.values())
        weights = np.array([weights[c] if c in weights else mean_weight for c in unique_states], dtype=float)
        weights/=weights.sum()
    else:
        weights = np.ones(nc, dtype=float)/nc

    # set up dummy matrix
    W = np.ones((nc,nc), dtype=float)

    mugration_GTR = GTR.custom(pi = weights, W=W, alphabet = np.array(alphabet))
    mugration_GTR.profile_map[missing_char] = np.ones(nc)
    mugration_GTR.ambiguous=missing_char

    ###########################################################################
    ### set up treeanc
    ###########################################################################
    treeanc = TreeAnc(params.tree, gtr=mugration_GTR, verbose=params.verbose)
    pseudo_seqs = [SeqRecord(id=n.name,name=n.name,
                   seq=Seq(reverse_alphabet[leaf_to_attr[n.name]] if n.name in leaf_to_attr else missing))
                   for n in treeanc.tree.get_terminals()]
    treeanc.aln = MultipleSeqAlignment(pseudo_seqs)

    treeanc.infer_ancestral_sequences(method='ml', infer_gtr=True,
            store_compressed=False, pc=5.0, marginal=True, normalized_rate=False,
            fixed_pi=weights if params.weights else None)


    ###########################################################################
    ### output
    ###########################################################################
    print("\nCompleted mugration model inference of attribute '%s' for"%attr,params.tree)

    bname = './'+os.path.basename(params.tree)
    gtr_name = bname + '.GTR.txt'
    with open(gtr_name, 'w') as ofile:
        ofile.write(str(treeanc.gtr)+'\n')
        print("\nSaved inferred mugration model as:", gtr_name)

    terminal_count = 0
    for n in treeanc.tree.find_clades():
        if n.up is None:
            continue
        n.confidence=None
        # due to a bug in older versions of biopython that truncated filenames in nexus export
        # we truncate them by hand and make them unique.
        if n.is_terminal() and len(n.name)>40 and bioversion<"1.69":
            n.name = n.name[:35]+'_%03d'%terminal_count
            terminal_count+=1
        n.comment= '&%s="'%attr + letter_to_state[n.sequence[0]] +'"'

    if params.confidence:
        conf_name = bname+'.confidence.csv'
        with open(conf_name, 'w') as ofile:
            ofile.write('#name, '+', '.join(unique_states)+'\n')
            for n in treeanc.tree.find_clades():
                ofile.write(n.name + ', '+', '.join([str(x) for x in n.marginal_profile[0]])+'\n')
        print("Saved table with ancestral state confidences as:", conf_name)

    # write tree to file
    outtree_name = bname+'.mugration.nexus'
    Phylo.write(treeanc.tree, outtree_name, 'nexus')
    print("Saved annotated tree as:",outtree_name)