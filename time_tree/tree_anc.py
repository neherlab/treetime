"""
Class defines simple tree object with basic interface methdos: reading and 
saving from/to files, initializing leaves with sequences from the alignment, 
making ancestral state inferrence
"""

from Bio import Phylo
from Bio import AlignIO
import numpy as np

class TreeAnc(object, Bio.Phylo.BaseTree):
    def __init__(self):
        pass
    
    # input stuff    
    @classmethod
    def from_file(cls, inf):
        pass
    
    @classmethod
    def _from_newick(cls, inf):
        pass
    
    @classmethod
    def _from_json(cls, inf):
        pass

    # ancestral state reconstruction 
    def set_seqs_to_leaves(self, aln):
        pass

    def reconstruct_anc(self, model):
        pass


