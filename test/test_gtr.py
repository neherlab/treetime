"""
Script to test functions of GTR class
"""

import unittest,os,sys
from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
from treetime import treeanc as ta

class TestGTR(unittest.TestCase):

    def _random_Q(self):
        self.gtr = ta.GTR(ta.alphabets['nuc'])


        Pi = 1.0*np.random.randint(0,100,size=(5))
        Pi /= Pi.sum()
        Pi = np.diagflat(Pi)

        W = 1.0*np.random.randint(0,100,size=(5,5)) # with gaps
        W = W+W.T
        np.fill_diagonal(W, 0)

        Wdiag = - ((W.T *  np.diagonal(Pi)).T).sum(0) / np.diagonal(Pi)
        np.fill_diagonal(W, Wdiag)

        Q = Pi.dot(W)
        assert((Q.sum(0) < 1e-10).sum() == 4)

        self.gtr.W = W
        self.gtr.Pi = Pi

        self.gtr.
    def test_reversibility(self):
        pass

    def test_jc(self):
        # can instantiate
        gtr = ta.GTR.standard('Jukes-Cantor')
        #concentrations are assigned correctly
        assert (gtr.Pi.sum() == 1.0)
        # the matrix is the rate matrix
        assert abs((gtr.Pi.dot(gtr.W)).sum(0).sum() < 1e-15)
        # eigendecomposition is made correctly
        a = gtr.v.shape(0)
        assert abs((gtr.v.dot(gtr.v_inv) - np.identity(a)).sum() < 1e-10)
        assert gtr.v.sum() > 1e-10 # **and** v is not zero

    def test_random(self):
        # can instantiate
        gtr = ta.GTR.standard('random')
        #concentrations are assigned correctly
        assert (gtr.Pi.sum() == 1.0)
        # the matrix is the rate matrix
        assert abs((gtr.Pi.dot(gtr.W)).sum(0).sum() < 1e-15)
        # eigendecomposition is made correctly
        a = gtr.v.shape(0)
        assert abs((gtr.v.dot(gtr.v_inv) - np.identity(a)).sum() < 1e-10)
        assert gtr.v.sum() > 1e-10 # **and** v is not zero


if __name__ == '__main__':





