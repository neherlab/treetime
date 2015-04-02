"""
Class, which contains methods to optimize branch lengths given the time
constraints set to leaves
"""
from tree_anc import TreeAnc
import numpy as np
from Bio import AlignIO


class TreeTime(TreeAnc):

    def __init__(self):
        pass

    def set_dates_to_leaves(self, aln):
        pass

    def dates_to_dist(self):
        pass

    def optimize_branch_len(self, model):
        pass
