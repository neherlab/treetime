from __future__ import print_function
import sys
if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO


# Tests
def import_short_test():
    print("testing short imports")
    from treetime import GTR
    from treetime import TreeTime
    from treetime import TreeAnc
    from treetime import seq_utils


def test_GTR():
    from treetime import GTR
    import numpy as np
    for model in ['Jukes-Cantor']:
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
    fasta = str(os.path.join(root_dir, 'treetime_examples/data/h3n2_na/h3n2_na_20.fasta'))
    nwk = str(os.path.join(root_dir, 'treetime_examples/data/h3n2_na/h3n2_na_20.nwk'))

    for marginal in [True, False]:
        print('loading flu example')
        t = TreeAnc(gtr='Jukes-Cantor', tree=nwk, aln=fasta)
        print('ancestral reconstruction' + ("marginal" if marginal else "joint"))
        t.reconstruct_anc(method='ml', marginal=marginal)
        assert "".join(t.tree.root.sequence) == 'ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATGCTTCTTCATGCAAATTGCCATCTTGATAACTACTGTAACATTGCATTTCAAGCAATATGAATTCAACTCCCCCCCAAACAACCAAGTGATGCTGTGTGAACCAACAATAATAGAAAGAAACATAACAGAGATAGTGTATCTGACCAACACCACCATAGAGAAGGAAATATGCCCCAAACCAGCAGAATACAGAAATTGGTCAAAACCGCAATGTGGCATTACAGGATTTGCACCTTTCTCTAAGGACAATTCGATTAGGCTTTCCGCTGGTGGGGACATCTGGGTGACAAGAGAACCTTATGTGTCATGCGATCCTGACAAGTGTTATCAATTTGCCCTTGGACAGGGAACAACACTAAACAACGTGCATTCAAATAACACAGTACGTGATAGGACCCCTTATCGGACTCTATTGATGAATGAGTTGGGTGTTCCTTTTCATCTGGGGACCAAGCAAGTGTGCATAGCATGGTCCAGCTCAAGTTGTCACGATGGAAAAGCATGGCTGCATGTTTGTATAACGGGGGATGATAAAAATGCAACTGCTAGCTTCATTTACAATGGGAGGCTTGTAGATAGTGTTGTTTCATGGTCCAAAGAAATTCTCAGGACCCAGGAGTCAGAATGCGTTTGTATCAATGGAACTTGTACAGTAGTAATGACTGATGGAAGTGCTTCAGGAAAAGCTGATACTAAAATACTATTCATTGAGGAGGGGAAAATCGTTCATACTAGCACATTGTCAGGAAGTGCTCAGCATGTCGAAGAGTGCTCTTGCTATCCTCGATATCCTGGTGTCAGATGTGTCTGCAGAGACAACTGGAAAGGCTCCAATCGGCCCATCGTAGATATAAACATAAAGGATCATAGCATTGTTTCCAGTTATGTGTGTTCAGGACTTGTTGGAGACACACCCAGAAAAAACGACAGCTCCAGCAGTAGCCATTGTTTGGATCCTAACAATGAAGAAGGTGGTCATGGAGTGAAAGGCTGGGCCTTTGATGATGGAAATGACGTGTGGATGGGAAGAACAATCAACGAGACGTCACGCTTAGGGTATGAAACCTTCAAAGTCATTGAAGGCTGGTCCAACCCTAAGTCCAAATTGCAGATAAATAGGCAAGTCATAGTTGACAGAGGTGATAGGTCCGGTTATTCTGGTATTTTCTCTGTTGAAGGCAAAAGCTGCATCAATCGGTGCTTTTATGTGGAGTTGATTAGGGGAAGAAAAGAGGAAACTGAAGTCTTGTGGACCTCAAACAGTATTGTTGTGTTTTGTGGCACCTCAGGTACATATGGAACAGGCTCATGGCCTGATGGGGCGGACCTCAATCTCATGCCTATA'

    print('testing LH normalization')
    from Bio import Phylo,AlignIO
    tiny_tree = Phylo.read(StringIO("((A:0.60100000009,B:0.3010000009):0.1,C:0.2):0.001;"), 'newick')
    tiny_aln = AlignIO.read(StringIO(">A\nAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT\n"
                                     ">B\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n"
                                     ">C\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"), 'fasta')

    mygtr = GTR.custom(alphabet = np.array(['A', 'C', 'G', 'T']), pi = np.array([0.9, 0.06, 0.02, 0.02]), W=np.ones((4,4)))
    t = TreeAnc(gtr=mygtr, tree=tiny_tree, aln=tiny_aln)
    t.reconstruct_anc('ml', marginal=True, debug=True)
    lhsum =  np.exp(t.sequence_LH(pos=np.arange(4**3))).sum()
    print (lhsum)
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
    import numpy as np
    try:
        from itertools import izip
    except ImportError:  #python3.x
        izip = zip
    from collections import defaultdict
    def exclusion(a, b):
        """
        Intersection of two lists
        """
        return list(set(a) - set(b))

    tiny_tree = Phylo.read(StringIO("((A:.060,B:.01200)C:.020,D:.0050)E:.004;"), 'newick')
    mygtr = GTR.custom(alphabet = np.array(['A', 'C', 'G', 'T']),
                       pi = np.array([0.15, 0.95, 0.05, 0.3]), W=np.ones((4,4)))
    seq = np.random.choice(mygtr.alphabet, p=mygtr.Pi, size=400)


    myTree = TreeAnc(gtr=mygtr, tree=tiny_tree, aln=None, verbose=4)

    # simulate evolution, set resulting sequence as ref_seq
    tree = myTree.tree
    seq_len = 400
    tree.root.ref_seq = np.random.choice(mygtr.alphabet, p=mygtr.Pi, size=seq_len)
    print ("Root sequence: " + ''.join(tree.root.ref_seq))
    mutation_list = defaultdict(list)
    for node in tree.find_clades():
        for c in node.clades:
            c.up = node
        if hasattr(node, 'ref_seq'):
            continue
        t = node.branch_length
        p = mygtr.evolve( seq_utils.seq2prof(node.up.ref_seq, mygtr.profile_map), t)
        # normalie profile
        p=(p.T/p.sum(axis=1)).T
        # sample mutations randomly
        ref_seq_idxs = np.array([int(np.random.choice(np.arange(p.shape[1]), p=p[k])) for k in np.arange(p.shape[0])])

        node.ref_seq = np.array([mygtr.alphabet[k] for k in ref_seq_idxs])

        node.ref_mutations = [(anc, pos, der) for pos, (anc, der) in
                            enumerate(izip(node.up.ref_seq, node.ref_seq)) if anc!=der]
        for anc, pos, der in node.ref_mutations:
            mutation_list[pos].append((node.name, anc, der))
        print (node.name, len(node.ref_mutations), node.ref_mutations)

    # set as the starting sequences to the terminal nodes:
    alnstr = ""
    i = 1
    for leaf in tree.get_terminals():
        alnstr += ">" + leaf.name + "\n" + ''.join(leaf.ref_seq) + '\n'
        i += 1
    print (alnstr)
    myTree.aln = AlignIO.read(StringIO(alnstr), 'fasta')
    myTree._attach_sequences_to_nodes()
    # reconstruct ancestral sequences:
    myTree._ml_anc_joint(debug=True)

    diff_count = 0
    mut_count = 0
    for node in myTree.tree.find_clades():
        if node.up is not None:
            mut_count += len(node.ref_mutations)
            diff_count += np.sum(node.sequence != node.ref_seq)==0
            if np.sum(node.sequence != node.ref_seq):
                print("%s: True sequence does not equal inferred sequence. parent %s"%(node.name, node.up.name))
            else:
                print("%s: True sequence equals inferred sequence. parent %s"%(node.name, node.up.name))
        print (node.name, np.sum(node.sequence != node.ref_seq), np.where(node.sequence != node.ref_seq), len(node.mutations), node.mutations)

    # the assignment of mutations to the root node is probabilistic. Hence some differences are expected
    assert diff_count/seq_len<2*(1.0*mut_count/seq_len)**2

    # prove the likelihood value calculation is correct
    LH = myTree.ancestral_likelihood()
    LH_p = (myTree.tree.sequence_LH)

    print ("Difference between reference and inferred LH:", (LH - LH_p).sum())
    assert ((LH - LH_p).sum())<1e-9

    return myTree


def test_seq_joint_lh_is_max():
    """
    For a single-char sequence, perform joint ancestral sequence reconstruction
    and prove that this reconstruction is the most likely one by comparing to all
    possible reconstruction variants (brute-force).
    """

    from treetime import TreeAnc, GTR
    from treetime import seq_utils
    from Bio import Phylo, AlignIO
    import numpy as np

    mygtr = GTR.custom(alphabet = np.array(['A', 'C', 'G', 'T']), pi = np.array([0.91, 0.05, 0.02, 0.02]), W=np.ones((4,4)))
    tiny_tree = Phylo.read(StringIO("((A:.0060,B:.30)C:.030,D:.020)E:.004;"), 'newick')

    #terminal node sequences (single nuc)
    A_char = 'A'
    B_char = 'C'
    D_char = 'G'

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

        myTree_1.reconstruct_anc(method='ml', marginal=False, debug=True)
        logLH = myTree_1.tree.sequence_LH
        return logLH

    ref = ref_lh()
    real  = real_lh()

    print(abs(ref.max() - real) )
    # joint chooses the most likely realization of the tree
    assert(abs(ref.max() - real) < 1e-10)
    return ref, real



