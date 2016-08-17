from __future__ import print_function


# Tests
def import_short_test():
    print("testing short imports")
    from treetime import GTR
    from treetime import TreeTime
    from treetime import TreeAnc


def test_GTR():
    from treetime import GTR
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


def test_ancestral():
    import os
    from Bio import AlignIO
    import numpy as np
    from treetime import TreeAnc, GTR
    root_dir = os.path.dirname(os.path.realpath(__file__))
    fasta = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.20.fasta')
    nwk = os.path.join(root_dir, '../data/H3N2_NA_allyears_NA.20.nwk')

    print('loading flu example')
    t = TreeAnc(gtr='Jukes-Cantor', tree=nwk, aln=fasta)

    print('ancestral reconstruction')
    t.reconstruct_anc(method='ml')

    assert "".join(t.tree.root.sequence) == 'ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATGCTTCTTCATGCAAATTGCCATCTTGATAACTACTGTAACATTGCATTTCAAGCAATATGAATTCAACTCCCCCCCAAACAACCAAGTGATGCTGTGTGAACCAACAATAATAGAAAGAAACATAACAGAGATAGTGTATCTGACCAACACCACCATAGAGAAGGAAATATGCCCCAAACCAGCAGAATACAGAAATTGGTCAAAACCGCAATGTGGCATTACAGGATTTGCACCTTTCTCTAAGGACAATTCGATTAGGCTTTCCGCTGGTGGGGACATCTGGGTGACAAGAGAACCTTATGTGTCATGCGATCCTGACAAGTGTTATCAATTTGCCCTTGGACAGGGAACAACACTAAACAACGTGCATTCAAATAACACAGTACGTGATAGGACCCCTTATCGGACTCTATTGATGAATGAGTTGGGTGTTCCTTTTCATCTGGGGACCAAGCAAGTGTGCATAGCATGGTCCAGCTCAAGTTGTCACGATGGAAAAGCATGGCTGCATGTTTGTATAACGGGGGATGATAAAAATGCAACTGCTAGCTTCATTTACAATGGGAGGCTTGTAGATAGTGTTGTTTCATGGTCCAAAGAAATTCTCAGGACCCAGGAGTCAGAATGCGTTTGTATCAATGGAACTTGTACAGTAGTAATGACTGATGGAAGTGCTTCAGGAAAAGCTGATACTAAAATACTATTCATTGAGGAGGGGAAAATCGTTCATACTAGCACATTGTCAGGAAGTGCTCAGCATGTCGAAGAGTGCTCTTGCTATCCTCGATATCCTGGTGTCAGATGTGTCTGCAGAGACAACTGGAAAGGCTCCAATCGGCCCATCGTAGATATAAACATAAAGGATCATAGCATTGTTTCCAGTTATGTGTGTTCAGGACTTGTTGGAGACACACCCAGAAAAAACGACAGCTCCAGCAGTAGCCATTGTTTGGATCCTAACAATGAAGAAGGTGGTCATGGAGTGAAAGGCTGGGCCTTTGATGATGGAAATGACGTGTGGATGGGAAGAACAATCAACGAGACGTCACGCTTAGGGTATGAAACCTTCAAAGTCATTGAAGGCTGGTCCAACCCTAAGTCCAAATTGCAGATAAATAGGCAAGTCATAGTTGACAGAGGTGATAGGTCCGGTTATTCTGGTATTTTCTCTGTTGAAGGCAAAAGCTGCATCAATCGGTGCTTTTATGTGGAGTTGATTAGGGGAAGAAAAGAGGAAACTGAAGTCTTGTGGACCTCAAACAGTATTGTTGTGTTTTGTGGCACCTCAGGTACATATGGAACAGGCTCATGGCCTGATGGGGCGGACCTCAATCTCATGCCTATA'

    print('testing LH normalization')
    from StringIO import StringIO
    from Bio import Phylo,AlignIO
    tiny_tree = Phylo.read(StringIO("((A:0.60100000009,B:0.3010000009):0.1,C:0.2):0.001;"), 'newick')
    tiny_aln = AlignIO.read(StringIO(">A\nAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT\n"
                                     ">B\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n"
                                     ">C\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"), 'fasta')

    mygtr = GTR.custom(alphabet = np.array(['A', 'C', 'G', 'T']), pi = np.array([0.9, 0.06, 0.02, 0.02]), W=np.ones((4,4)))
    t = TreeAnc(gtr=mygtr, tree=tiny_tree, aln=tiny_aln)
    t.reconstruct_anc('ml', marginal=True)
    lhsum = (np.exp(t.tree.root.lh_prefactor)*t.tree.root.profile.sum(axis=1)).sum()
    assert(np.abs(lhsum-1.0)<1e-6)

    t.optimize_branch_len()