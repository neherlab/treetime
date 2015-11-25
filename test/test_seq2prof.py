"""
Script to test functions of TreeAnc module
"""

import unittest,os,sys
from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from treetime import treeanc as ta

class TestSeq2Prof(unittest.TestCase):

    def _gen_random_seq(self, L, tp='nuc'):
        s = np.random.choice(ta._full_nc_profile.keys(), L)
        return s

    def test_seq2prof(self,):
        s = self._gen_random_seq(1000)
        p = ta.seq2prof(s, 'nuc')
        p_ref = np.array([ta._full_nc_profile[k] for k in s])
        assert (p_ref - p).sum() < 1e-16


