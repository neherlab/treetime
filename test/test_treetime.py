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
        print('Frequency sum:', myGTR.Pi.sum())
        assert (myGTR.Pi.sum() - 1.0)**2<1e-15
        # the matrix is the rate matrix
        assert abs((myGTR.Pi.dot(myGTR.W)).sum(0).sum() < 1e-15)
        # eigendecomposition is made correctly
        n_states = myGTR.v.shape[0]
        assert abs((myGTR.v.dot(myGTR.v_inv) - np.identity(n_states)).sum() < 1e-10)
        assert np.abs(myGTR.v.sum()) > 1e-10 # **and** v is not zero


def test_TreeTime():
    import os
    from Bio import AlignIO
    import numpy as np
    root_dir = os.path.dirname(os.path.realpath(__file__))
    fasta = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.20.fasta')
    nwk = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.20.nwk')
    mdf = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.20.metadata.csv')

    print('loading flu example')
    from treetime.treetime.gtr import GTR
    from treetime.treetime import io
    gtr = GTR.standard()
    t = io.treetime_from_newick(gtr, nwk)

    print('attaching leaves')
    io.set_seqs_to_leaves(t, AlignIO.read(fasta, 'fasta'))
    io.read_metadata(t, mdf)
    print('rerooting and infer gtr')
    t.reroot_to_best_root(infer_gtr=True)
    print('time tree inference')
    #t.ml_t()

    assert "".join(t.tree.root.sequence) == 'ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTGCCACAATATGCTTCCTTATGCAAATTGCCATCCTGGTAACTACTGTAACATTGCATTTCAAGCAATATGAATGCAACTCCCCCCCAAACAACCAAGTGATGCTGTGTGAACCAACAATAATAGAAAGAAACATAACAGAGATAGTGTATCTGACCAACACCACCATAGAGAAGGAAATATGCCCCAAACTAGCAGAATACAGAAATTGGTCAAAGCCGCAATGTAACATTACAGGATTTGCACCTTTTTCTAAGGACAATTCGATTCGGCTTTCCGCTGGTGGGGACATCTGGGTGACAAGAGAACCTTATGTGTCATGCGATCCTGACAAGTGTTATCAATTTGCCCTTGGACAGGGAACAACACTAAACAACGGGCATTCAAATGACACAGTACATGATAGGACCCCTTATCGGACCCTATTGATGAATGAGTTGGGTGTTCCATTTCATTTGGGAACCAAGCAAGTGTGCATAGCATGGTCCAGCTCAAGTTGTCACGATGGAAAAGCATGGCTGCATGTTTGTGTAACGGGGGATGATGAAAATGCAACTGCTAGCTTCATTTACAATGGGAGGCTTGTAGATAGTATTGGTTCATGGTCCAAAAAAATCCTCAGGACCCAGGAGTCGGAATGCGTTTGTATCAATGGAACTTGTACAGTGGTAATGACTGATGGGAGTGCTTCAGGAAAAGCTGATACTAAAATACTATTCATTGAGGAGGGGAAAATCGTTCATACTAGCACATTGTCAGGAAGTGCTCAGCATGTCGAGGAGTGCTCCTGTTATCCTCGATATCCTGGTGTCAGATGTGTCTGCAGAGACAACTGGAAAGGCTCCAATAGGCCCATCGTAGATATAAATGTAAAGGATTATAGCATTGTTTCCAGTTATGTGTGCTCAGGACTTGTTGGAGACACACCCAGAAAAAACGACAGCTCCAGCAGTAGCCATTGCTTGGATCCTAACAATGAGGAAGGTGGTCATGGAGTGAAAGGCTGGGCCTTTGATGATGGAAATGACGTGTGGATGGGAAGAACGATCAGCGAGAAGTTACGCTCAGGATATGAAACCTTCAAAGTCATTGAAGGCTGGTCCAAACCTAACTCCAAATTGCAGATAAATAGGCAAGTCATAGTTGACAGAGGTAATAGGTCCGGTTATTCTGGTATTTTCTCTGTTGAAGGCAAAAGCTGCATCAATCGGTGCTTTTATGTGGAGTTGATAAGGGGAAGGAAACAGGAAACTGAAGTCTTGTGGACCTCAAACAGTATTGTTGTGTTTTGTGGCACCTCAGGTACATATGGAACAGGCTCATGGCCTGATGGGGCGGACATCAATCTCATGCCTATA'
    #assert np.abs(t.tree.root.numdate - 1997)<1
