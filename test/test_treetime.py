from __future__ import print_function


# Tests
def import_test():
    print("testing imports")
    from treetime.treetime.gtr import GTR
    from treetime.treetime.treetime import TreeTime
    from treetime.treetime.treeanc import TreeAnc
    from treetime.treetime import io, utils

# Tests
def import_test():
    print("testing short imports")
    from treetime import GTR
    from treetime import TreeTime
    from treetime import TreeAnc


def test_GTR():
    from treetime.treetime.gtr import GTR
    import numpy as np
    for model in ['Jukes-Cantor', 'random']:
        print('testing GTR, model:',model)
        myGTR = GTR.standard(model, alphabet='nuc')
        assert (myGTR.Pi.sum() == 1.0)
        # the matrix is the rate matrix
        assert abs((myGTR.Pi.dot(myGTR.W)).sum(0).sum() < 1e-15)
        # eigendecomposition is made correctly
        n_states = myGTR.v.shape[0]
        assert abs((myGTR.v.dot(myGTR.v_inv) - np.identity(n_states)).sum() < 1e-10)
        assert np.abs(myGTR.v.sum()) > 1e-10 # **and** v is not zero
