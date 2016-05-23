from __future__ import print_function, division
import numpy as np
import os,sys,copy,json
from Bio import Phylo, AlignIO
import treetime

#  set the configuration parameters
reuse_branch_len=True  # Do we trust the branch lengths of the tree given? 
infer_gtr=False  # Infer GTR model from the data given (Use default if False)?
reroot = False  # Find the root node which maximizes the molecular clock correlation?

use_mu=True  # Should we use the mutation rate (if False, it will be inferred from the molecular clock)? 
mu=1e-2  #  #mutations per years per position
resolve_poly=True  # Try to resolve multiple mergers
coalescent = True  # Run the coalescent model (takes care about polytomies resolution)
Tc = 0.01  #  coalescent distance (%difference )
relax_mu = True  #  Relax mutation rate? 
slack = 0.0  # How much we allow for variation between parents and children (from 0 to 1)
coupling = 0.0  # How much we allow for variation between sister branches (from 0 to 1) 

#root = os.path.join(__dirname__ , "../data/") #  where to get data from
root = "../data" #= os.path.join(__dirname__ , "../data/") #  where to get data from

# specify filenames
nwk = os.path.join(root, "H3N2_NA_allyears_NA.20.nwk")
aln = os.path.join(root, "H3N2_NA_allyears_NA.20.fasta")
meta = os.path.join(root, "H3N2_NA_allyears_NA.20.metadata.csv")

if __name__=="__main__":


    gtr = treetime.GTR.standard() # always start from standard J-C model
    tt = treetime.treetime_from_newick(gtr, nwk)  # load tree
    treetime.set_seqs_to_leaves(tt, AlignIO.read(aln, 'fasta'))  # assign sequences
    treetime.read_metadata(tt, meta)  # load meta data 

    ##  Pipeline
    if not reuse_branch_len:
        tt.optimize_seq_and_branch_len(False, True)
    
    if reroot:
        # This function does the following:
        # 1. re-root 
        # 2. infer GTR model if needed
        # 3. ancestral inference
        # 4. init time constraints
        tt.reroot_to_best_root(infer_gtr=infer_gtr) # we can 

    elif infer_gtr: # 
        tt.infer_gtr()  # will make first shot in ancestral state reconstruction
        tt.init_date_constraints(slope=mu, ancestral_inference=True) 
    else:
        # just initialize  the objects needed to run the TreeTime-ML, 
        # and infer ancestral sequences
        tt.init_date_constraints(slope=mu, ancestral_inference=True)


    # run TreeTime-ML optimization
    tt.ml_t()

    # some post-processing
    if coalescent:
        tt.coalescent_model(Tc) # Run the coalescent model + resolve polytomies
    elif resolve_poly:
        tt.resolve_polytomies() # Resolve multiple mergers
    else:
        pass  # do nothing

    if relax_mu:  # Relaxed molecular clock 
        tt.relaxed_clock(slack, coupling)



