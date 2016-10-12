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
        assert (myGTR.Pi.sum() - 1.0)**2<1e-14
        # the matrix is the rate matrix
        assert abs(myGTR.Q.sum(0)).sum() < 1e-14
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


def test_seq_joint_reconstruction_correct():
    """
    evolve the random sequence, get the alignment at the leaf nodes. 
    Reconstruct the sequences of the internal nodes (joint)
    and prove the reconstruction is correct.
    In addition, compute the likelihood of the particular realization of the 
    sequences on the tree and prove that this likelihood is exactly the same 
    as calculated in the joint reconstruction
    """

    from treetime import TreeAnc, GTR
    from treetime import seq_utils
    from Bio import Phylo, AlignIO
    from StringIO import StringIO
    import numpy as np
    tiny_tree = Phylo.read(StringIO("((A:.060,B:.30)C:.30,D:.020)E:.04;"), 'newick')
    mygtr = GTR.custom(alphabet = np.array(['A', 'C', 'G', 'T']), pi = np.array([0.15, 0.95, 0.05, 0.3]), W=np.ones((4,4)))
    seq = np.random.choice(mygtr.alphabet, p=mygtr.Pi, size=400)


    myTree = TreeAnc(gtr=mygtr, tree=tiny_tree, aln=None, verbose=4)

    # simulate evolution, set resulting sequence as ref_seq
    tree = myTree.tree 
    tree.root.ref_seq = np.random.choice(mygtr.alphabet, p=mygtr.Pi, size=400)
    print ("Root sequence: " + ''.join(tree.root.ref_seq))
    for node in tree.find_clades():
        for c in node.clades: 
            c.up = node
        if hasattr(node, 'ref_seq'):
            continue
        t = node.branch_length
        p = mygtr.propagate_profile( seq_utils.seq2prof(node.up.ref_seq, mygtr.profile_map), t)
        # normalie profile
        p=(p.T/p.sum(axis=1)).T
        # smaple mutations randomly
        ref_seq_idxs = np.array([int(np.random.choice(np.arange(p.shape[1]), p=p[k])) for k in np.arange(p.shape[0])])
        # sample by most-likely
        ref_seq_idxs_ref = p.argmax(axis = 1)
        
        ref_seq = np.array([mygtr.alphabet[k] for k in ref_seq_idxs_ref])
        node.ref_seq = ref_seq

    # set as the starting sequencees to the terminal nodes:
    alnstr = ""
    i = 1
    for leaf in tree.get_terminals():
        alnstr += ">" + leaf.name + "\n" + ''.join(leaf.ref_seq) + '\n'
        i += 1 
    print (alnstr)
    myTree.aln = AlignIO.read(StringIO(alnstr), 'fasta')
    myTree.attach_sequences_to_nodes()
    # reconstruct ancestral sequences:
    myTree._ml_anc_joint( sample_from_profile='root')

    for node in myTree.tree.get_nonterminals():
        print (node.name + ":  " + str(np.sum(node.sequence != node.ref_seq)))

    # prove the likelihood valu calculation is correct
    LH = myTree.ancestral_likelihood()
    assert(abs(LH.sum() - myTree.tree.joint_seq_LH) < 1e-9)

    return myTree


def test_seq_joint_lh_is_max():
    """
    For a single-cahr sequence, make the joint ancestracl sequence reconstruction
    and prove that this reconstruction is the most likely one by comparing to all 
    possible reconstruction variants (brute-force). 
    """

    from treetime import TreeAnc, GTR
    from treetime import seq_utils
    from Bio import Phylo, AlignIO
    from StringIO import StringIO
    import numpy as np

    mygtr = GTR.custom(alphabet = np.array(['A', 'C', 'G', 'T']), pi = np.array([0.91, 0.05, 0.02, 0.02]), W=np.ones((4,4)))
    tiny_tree = Phylo.read(StringIO("((A:.0060,B:.30)C:.030,D:.020)E:.004;"), 'newick')
    
    #terminal node sequences (single nuc)    
    A_char = 'A'
    B_char = 'C'
    D_char = 'C'

    # for brute-force, expand them to the strings
    A_seq = ''.join(np.repeat(A_char,16))
    B_seq = ''.join(np.repeat(B_char,16))
    D_seq = ''.join(np.repeat(D_char,16))

    # 
    def ref_lh():
        """
        reference likelihood - LH values for all possible variants 
        of the internal node sequences
        """

        tiny_aln = AlignIO.read(StringIO(">A\n" + A_seq + "\n"
                                         ">B\n" + B_seq + "\n"
                                         ">D\n" + D_seq + "\n"
                                         ">C\nAAAACCCCGGGGTTTT\n"
                                         ">E\nACGTACGTACGTACGT\n"), 'fasta')
    
        myTree = TreeAnc(gtr=mygtr, tree = tiny_tree,
                            aln =tiny_aln, verbose = 4)
    
        logLH_ref = myTree.ancestral_likelihood()

        return logLH_ref

    # 
    def real_lh():
        """
        Likelihood of the sequences calculated by the joint ancestral 
        sequence reconstruction
        """
        tiny_aln_1 = AlignIO.read(StringIO(">A\n"+A_char+"\n"
                                           ">B\n"+B_char+"\n"
                                           ">D\n"+D_char+"\n"), 'fasta')
        
        myTree_1 = TreeAnc(gtr=mygtr, tree = tiny_tree,
                            aln=tiny_aln_1, verbose = 4)
        
        myTree_1._ml_anc_joint()
        logLH = myTree_1.tree.joint_seq_LH
        return logLH
    
    ref = ref_lh()
    real  = real_lh()

    # joint chooses the most likely realization of the tree
    assert(abs(ref.max() - real) < 1e-10)
    return ref, real



