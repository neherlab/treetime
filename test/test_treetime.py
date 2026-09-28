from __future__ import print_function
from io import StringIO
import pytest


# Tests
def test_import_short():
    print('testing short imports')
    from treetime import GTR
    from treetime import TreeTime
    from treetime import TreeAnc
    from treetime import seq_utils


def test_assign_gamma(root_dir=None):
    print('testing assign gamma')
    import os
    from treetime import TreeTime
    from treetime.utils import parse_dates

    if root_dir is None:
        root_dir = os.path.dirname(os.path.realpath(__file__))
    fasta = root_dir + '/treetime_examples/data/h3n2_na/h3n2_na_20.fasta'
    nwk = root_dir + '/treetime_examples/data/h3n2_na/h3n2_na_20.nwk'
    dates = parse_dates(root_dir + '/treetime_examples/data/h3n2_na/h3n2_na_20.metadata.csv')
    seq_kwargs = {'marginal_sequences': True, 'branch_length_mode': 'input', 'reconstruct_tip_states': False}
    tt_kwargs = {'clock_rate': 0.0001, 'time_marginal': 'assign'}
    myTree = TreeTime(
        gtr='Jukes-Cantor',
        tree=nwk,
        use_fft=False,
        aln=fasta,
        verbose=1,
        dates=dates,
        precision=3,
        debug=True,
        rng_seed=1234,
    )

    def assign_gamma(tree):
        return tree

    success = myTree.run(infer_gtr=False, assign_gamma=assign_gamma, max_iter=1, verbose=3, **seq_kwargs, **tt_kwargs)
    assert success


def test_GTR(root_dir=None):
    from treetime import GTR
    import numpy as np
    import os

    if root_dir is None:
        root_dir = os.path.dirname(os.path.realpath(__file__))
    ##check custom GTR model
    custom_gtr = root_dir + '/test_sequence_evolution_model.txt'
    gtr = GTR.from_file(custom_gtr)
    assert (gtr.Pi.sum() - 1.0) ** 2 < 1e-14
    assert np.allclose(gtr.Pi, np.array([0.3088, 0.1897, 0.2335, 0.2581, 0.0099]))
    assert np.all(gtr.alphabet == np.array(['A', 'C', 'G', 'T', '-']))
    assert abs(gtr.mu - 1.0) < 1e-4
    assert abs(gtr.Q.sum(0)).sum() < 1e-14
    assert np.allclose(
        gtr.W,
        np.array(
            [
                [0, 0.7003, 3.0669, 0.2651, 0.9742],
                [0.7003, 0, 0.3354, 3.399, 0.999],
                [3.0669, 0.3354, 0, 0.4258, 0.9892],
                [0.2651, 3.399, 0.4258, 0, 0.9848],
                [0.9742, 0.999, 0.9892, 0.9848, 0],
            ]
        ),
        atol=1e-4,
    )
    for model in ['Jukes-Cantor']:
        print('testing GTR, model:', model)
        myGTR = GTR.standard(model, alphabet='nuc')
        print('Frequency sum:', myGTR.Pi.sum())
        assert (myGTR.Pi.sum() - 1.0) ** 2 < 1e-14
        # the matrix is the rate matrix
        assert abs(myGTR.Q.sum(0)).sum() < 1e-14
        # eigendecomposition is made correctly
        n_states = myGTR.v.shape[0]
        assert abs((myGTR.v.dot(myGTR.v_inv) - np.identity(n_states)).sum() < 1e-10)
        assert np.abs(myGTR.v.sum()) > 1e-10  # **and** v is not zero


def test_reconstruct_discrete_traits():
    from Bio import Phylo
    from treetime.wrappers import reconstruct_discrete_traits

    # Create a minimal tree with traits to reconstruct.
    tiny_tree = Phylo.read(StringIO('((A:0.60100000009,B:0.3010000009):0.1,C:0.2):0.001;'), 'newick')
    traits = {
        'A': '?',
        'B': 'North America',
        'C': 'West Asia',
    }

    # Reconstruct traits with "?" as missing data.
    mugration, letter_to_state, reverse_alphabet = reconstruct_discrete_traits(
        tiny_tree,
        traits,
        missing_data='?',
    )

    # With two known states, the letters "A" and "B" should be in the alphabet
    # mapping to those states.
    assert 'A' in letter_to_state
    assert 'B' in letter_to_state

    # The letter for missing data should be the next letter in the alphabet,
    # following the two known state letters.
    assert letter_to_state['C'] == '?'


def test_ancestral(root_dir=None):
    import os
    from Bio import AlignIO
    import numpy as np
    from treetime import TreeAnc, GTR

    if root_dir is None:
        root_dir = os.path.dirname(os.path.realpath(__file__))
    fasta = str(os.path.join(root_dir, 'treetime_examples/data/h3n2_na/h3n2_na_20.fasta'))
    nwk = str(os.path.join(root_dir, 'treetime_examples/data/h3n2_na/h3n2_na_20.nwk'))

    for marginal in [True, False]:
        print('loading flu example')
        t = TreeAnc(gtr='Jukes-Cantor', tree=nwk, aln=fasta, rng_seed=1234)
        print('ancestral reconstruction' + ('marginal' if marginal else 'joint'))
        t.reconstruct_anc(method='ml', marginal=marginal)
        assert (
            t.data.compressed_to_full_sequence(t.tree.root.cseq, as_string=True)
            == 'ATGAATCCAAATCAAAAGATAATAACGATTGGCTCTGTTTCTCTCACCATTTCCACAATATGCTTCTTCATGCAAATTGCCATCTTGATAACTACTGTAACATTGCATTTCAAGCAATATGAATTCAACTCCCCCCCAAACAACCAAGTGATGCTGTGTGAACCAACAATAATAGAAAGAAACATAACAGAGATAGTGTATCTGACCAACACCACCATAGAGAAGGAAATATGCCCCAAACCAGCAGAATACAGAAATTGGTCAAAACCGCAATGTGGCATTACAGGATTTGCACCTTTCTCTAAGGACAATTCGATTAGGCTTTCCGCTGGTGGGGACATCTGGGTGACAAGAGAACCTTATGTGTCATGCGATCCTGACAAGTGTTATCAATTTGCCCTTGGACAGGGAACAACACTAAACAACGTGCATTCAAATAACACAGTACGTGATAGGACCCCTTATCGGACTCTATTGATGAATGAGTTGGGTGTTCCTTTTCATCTGGGGACCAAGCAAGTGTGCATAGCATGGTCCAGCTCAAGTTGTCACGATGGAAAAGCATGGCTGCATGTTTGTATAACGGGGGATGATAAAAATGCAACTGCTAGCTTCATTTACAATGGGAGGCTTGTAGATAGTGTTGTTTCATGGTCCAAAGAAATTCTCAGGACCCAGGAGTCAGAATGCGTTTGTATCAATGGAACTTGTACAGTAGTAATGACTGATGGAAGTGCTTCAGGAAAAGCTGATACTAAAATACTATTCATTGAGGAGGGGAAAATCGTTCATACTAGCACATTGTCAGGAAGTGCTCAGCATGTCGAAGAGTGCTCTTGCTATCCTCGATATCCTGGTGTCAGATGTGTCTGCAGAGACAACTGGAAAGGCTCCAATCGGCCCATCGTAGATATAAACATAAAGGATCATAGCATTGTTTCCAGTTATGTGTGTTCAGGACTTGTTGGAGACACACCCAGAAAAAACGACAGCTCCAGCAGTAGCCATTGTTTGGATCCTAACAATGAAGAAGGTGGTCATGGAGTGAAAGGCTGGGCCTTTGATGATGGAAATGACGTGTGGATGGGAAGAACAATCAACGAGACGTCACGCTTAGGGTATGAAACCTTCAAAGTCATTGAAGGCTGGTCCAACCCTAAGTCCAAATTGCAGATAAATAGGCAAGTCATAGTTGACAGAGGTGATAGGTCCGGTTATTCTGGTATTTTCTCTGTTGAAGGCAAAAGCTGCATCAATCGGTGCTTTTATGTGGAGTTGATTAGGGGAAGAAAAGAGGAAACTGAAGTCTTGTGGACCTCAAACAGTATTGTTGTGTTTTGTGGCACCTCAGGTACATATGGAACAGGCTCATGGCCTGATGGGGCGGACCTCAATCTCATGCCTATA'
        )

    print('testing LH normalization')
    from Bio import Phylo, AlignIO

    tiny_tree = Phylo.read(StringIO('((A:0.60100000009,B:0.3010000009):0.1,C:0.2):0.001;'), 'newick')
    tiny_aln = AlignIO.read(
        StringIO(
            '>A\nAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT\n'
            '>B\nAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT\n'
            '>C\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n'
        ),
        'fasta',
    )

    mygtr = GTR.custom(alphabet=np.array(['A', 'C', 'G', 'T']), pi=np.array([0.9, 0.06, 0.02, 0.02]), W=np.ones((4, 4)))
    t = TreeAnc(gtr=mygtr, tree=tiny_tree, aln=tiny_aln, rng_seed=1234)
    t.reconstruct_anc('ml', marginal=True, debug=True)
    lhsum = np.exp(t.sequence_LH(pos=np.arange(4**3))).sum()
    print(lhsum)
    assert np.abs(lhsum - 1.0) < 1e-6

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
    from collections import defaultdict

    tiny_tree = Phylo.read(StringIO('((A:.060,B:.01200)C:.020,D:.0050)E:.004;'), 'newick')
    mygtr = GTR.custom(alphabet=np.array(['A', 'C', 'G', 'T']), pi=np.array([0.15, 0.95, 0.05, 0.3]), W=np.ones((4, 4)))

    myTree = TreeAnc(gtr=mygtr, tree=tiny_tree, aln=None, verbose=4, rng_seed=1234)

    # simulate evolution, set resulting sequence as ref_seq
    tree = myTree.tree
    seq_len = 400
    tree.root.ref_seq = myTree.rng.choice(mygtr.alphabet, p=mygtr.Pi, size=seq_len)
    print('Root sequence: ' + ''.join(tree.root.ref_seq.astype('U')))
    mutation_list = defaultdict(list)
    for node in tree.find_clades():
        for c in node.clades:
            c.up = node
        if hasattr(node, 'ref_seq'):
            continue
        t = node.branch_length
        p = mygtr.evolve(seq_utils.seq2prof(node.up.ref_seq, mygtr.profile_map), t)
        # normalize profile
        p = (p.T / p.sum(axis=1)).T
        # sample mutations randomly
        ref_seq_idxs = np.array([int(myTree.rng.choice(np.arange(p.shape[1]), p=p[k])) for k in np.arange(p.shape[0])])

        node.ref_seq = np.array([mygtr.alphabet[k] for k in ref_seq_idxs])

        node.ref_mutations = [
            (anc, pos, der) for pos, (anc, der) in enumerate(zip(node.up.ref_seq, node.ref_seq)) if anc != der
        ]
        for anc, pos, der in node.ref_mutations:
            mutation_list[pos].append((node.name, anc, der))
        print(node.name, len(node.ref_mutations), node.ref_mutations)

    # set as the starting sequences to the terminal nodes:
    alnstr = ''
    i = 1
    for leaf in tree.get_terminals():
        alnstr += '>' + leaf.name + '\n' + ''.join(leaf.ref_seq.astype('U')) + '\n'
        i += 1
    print(alnstr)
    myTree.aln = AlignIO.read(StringIO(alnstr), 'fasta')
    # reconstruct ancestral sequences:
    myTree.infer_ancestral_sequences(final=True, debug=True, reconstruct_leaves=True)

    diff_count = 0
    mut_count = 0
    for node in myTree.tree.find_clades():
        if node.up is not None:
            mut_count += len(node.ref_mutations)
            diff_count += np.sum(node.sequence != node.ref_seq)
            if np.sum(node.sequence != node.ref_seq):
                print('%s: True sequence does not equal inferred sequence. parent %s' % (node.name, node.up.name))
            else:
                print('%s: True sequence equals inferred sequence. parent %s' % (node.name, node.up.name))

    # the assignment of mutations to the root node is probabilistic. Hence some differences are expected
    assert diff_count / seq_len < 2 * (1.0 * mut_count / seq_len) ** 2

    # prove the likelihood value calculation is correct
    LH = myTree.ancestral_likelihood()
    LH_p = myTree.tree.sequence_LH

    print('Difference between reference and inferred LH:', (LH - LH_p).sum())
    assert ((LH - LH_p).sum()) < 1e-9


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

    mygtr = GTR.custom(
        alphabet=np.array(['A', 'C', 'G', 'T']), pi=np.array([0.91, 0.05, 0.02, 0.02]), W=np.ones((4, 4))
    )
    tiny_tree = Phylo.read(StringIO('((A:.0060,B:.30)C:.030,D:.020)E:.004;'), 'newick')

    # terminal node sequences (single nuc)
    A_char = 'A'
    B_char = 'C'
    D_char = 'G'

    # for brute-force, expand them to the strings
    A_seq = ''.join(np.repeat(A_char, 16))
    B_seq = ''.join(np.repeat(B_char, 16))
    D_seq = ''.join(np.repeat(D_char, 16))

    #
    def ref_lh():
        """
        reference likelihood - LH values for all possible variants
        of the internal node sequences
        """

        tiny_aln = AlignIO.read(
            StringIO(
                '>A\n' + A_seq + '\n>B\n' + B_seq + '\n>D\n' + D_seq + '\n>C\nAAAACCCCGGGGTTTT\n>E\nACGTACGTACGTACGT\n'
            ),
            'fasta',
        )

        myTree = TreeAnc(gtr=mygtr, tree=tiny_tree, aln=tiny_aln, verbose=4, rng_seed=1234)

        logLH_ref = myTree.ancestral_likelihood()

        return logLH_ref

    #
    def real_lh():
        """
        Likelihood of the sequences calculated by the joint ancestral
        sequence reconstruction
        """
        tiny_aln_1 = AlignIO.read(StringIO('>A\n' + A_char + '\n>B\n' + B_char + '\n>D\n' + D_char + '\n'), 'fasta')

        myTree_1 = TreeAnc(gtr=mygtr, tree=tiny_tree, aln=tiny_aln_1, verbose=4, rng_seed=1234)

        myTree_1.reconstruct_anc(method='ml', marginal=False, debug=True)
        logLH = myTree_1.tree.sequence_LH
        return logLH

    ref = ref_lh()
    real = real_lh()

    print(abs(ref.max() - real))
    # joint chooses the most likely realization of the tree
    assert abs(ref.max() - real) < 1e-10


def test_marginal_mode_used_in_all_iterations(root_dir=None):
    """Verify that marginal reconstruction is used in every call to
    infer_ancestral_sequences when branch_length_mode='marginal'.

    Regression test for https://github.com/neherlab/treetime/issues/601:
    later iterations passed marginal_sequences via **seq_kwargs but not as the
    explicit `marginal` parameter, causing infer_ancestral_sequences to default
    to joint reconstruction.
    """
    import os
    from unittest.mock import patch
    from treetime import TreeTime
    from treetime.utils import parse_dates

    if root_dir is None:
        root_dir = os.path.dirname(os.path.realpath(__file__))

    fasta = root_dir + '/treetime_examples/data/h3n2_na/h3n2_na_20.fasta'
    nwk = root_dir + '/treetime_examples/data/h3n2_na/h3n2_na_20.nwk'
    dates = parse_dates(root_dir + '/treetime_examples/data/h3n2_na/h3n2_na_20.metadata.csv')

    myTree = TreeTime(
        gtr='Jukes-Cantor',
        tree=nwk,
        use_fft=False,
        aln=fasta,
        verbose=1,
        dates=dates,
        precision=3,
        debug=True,
        rng_seed=1234,
    )

    marginal_args_log = []
    original = myTree.infer_ancestral_sequences.__func__

    def tracking_wrapper(self, *args, marginal=False, **kwargs):
        marginal_args_log.append(marginal)
        return original(self, *args, marginal=marginal, **kwargs)

    with patch.object(type(myTree), 'infer_ancestral_sequences', tracking_wrapper):
        myTree.run(
            branch_length_mode='marginal',
            infer_gtr=False,
            max_iter=2,
            time_marginal=False,
        )

    assert len(marginal_args_log) > 0, "infer_ancestral_sequences was never called"
    for i, was_marginal in enumerate(marginal_args_log):
        assert was_marginal, (
            f"Call #{i} to infer_ancestral_sequences used marginal=False (joint) "
            f"when marginal reconstruction was requested"
        )


# The root ROOT has one child, MRCA, the common ancestor of all samples.
# Regression data for https://github.com/neherlab/treetime/issues/959
SINGLE_CHILD_ROOT_TREE = (
    '((((A:0.010,B:0.012)AB:0.004,(C:0.015,D:0.013)CD:0.003)ABCD:0.002,'
    '((E:0.011,F:0.016)EF:0.002,(G:0.014,H:0.009)GH:0.005)EFGH:0.003)MRCA:0.001)ROOT:0;'
)
SINGLE_CHILD_ROOT_DATES = {
    'A': 2014.0,
    'B': 2016.0,
    'C': 2019.0,
    'D': 2017.0,
    'E': 2013.0,
    'F': 2018.0,
    'G': 2019.0,
    'H': 2014.0,
}


def _single_child_root_tree_time(dates):
    from Bio import Phylo
    from treetime import TreeTime

    return TreeTime(
        tree=Phylo.read(StringIO(SINGLE_CHILD_ROOT_TREE), 'newick'),
        dates=dates,
        seq_len=1000,
        gtr='Jukes-Cantor',
        verbose=0,
        rng_seed=1234,
    )


def _write_single_child_root_inputs(tmp_path, dates):
    tree_file = tmp_path / 'tree.nwk'
    tree_file.write_text(SINGLE_CHILD_ROOT_TREE + '\n')
    dates_file = tmp_path / 'dates.tsv'
    dates_file.write_text('name\tdate\n' + ''.join(f'{name}\t{date}\n' for name, date in dates.items()))
    return tree_file, dates_file


def _run_cli(args):
    import matplotlib

    matplotlib.use('Agg')
    from treetime import make_parser

    params = make_parser().parse_args([str(arg) for arg in args])
    return params.func(params)


@pytest.mark.parametrize('root', ['least-squares', 'min_dev', 'oldest', 'A'])
def test_reroot_removes_undated_single_child_root(root):
    """Rerooting must not turn an undated single-child root into an extra leaf."""
    tt = _single_child_root_tree_time(SINGLE_CHILD_ROOT_DATES)

    tt.reroot(root)

    assert sorted(n.name for n in tt.tree.get_terminals()) == sorted(SINGLE_CHILD_ROOT_DATES)
    assert 'ROOT' not in {n.name for n in tt.tree.find_clades()}


def test_reroot_keeps_dated_single_child_root():
    """A dated single-child root carries data: after rerooting it is a dated leaf."""
    tt = _single_child_root_tree_time({**SINGLE_CHILD_ROOT_DATES, 'ROOT': 2008.0})

    tt.reroot('least-squares')

    leaves = {n.name: n for n in tt.tree.get_terminals()}
    assert sorted(leaves) == sorted([*SINGLE_CHILD_ROOT_DATES, 'ROOT'])
    assert leaves['ROOT'].raw_date_constraint == 2008.0


def test_reroot_to_current_single_child_root_keeps_tree():
    tt = _single_child_root_tree_time(SINGLE_CHILD_ROOT_DATES)

    tt.reroot(tt.tree.root)

    assert tt.tree.root.name == 'ROOT'
    assert [c.name for c in tt.tree.root] == ['MRCA']
    assert sorted(n.name for n in tt.tree.get_terminals()) == sorted(SINGLE_CHILD_ROOT_DATES)


def test_plot_root_to_tip_after_reroot_of_single_child_root():
    import matplotlib

    matplotlib.use('Agg')
    tt = _single_child_root_tree_time(SINGLE_CHILD_ROOT_DATES)
    tt.reroot('least-squares')

    tt.plot_root_to_tip()


@pytest.mark.parametrize(
    ('command', 'output_tree', 'output_format'),
    [([], 'timetree.nexus', 'nexus'), (['clock'], 'rerooted.newick', 'newick')],
)
def test_cli_reroot_of_single_child_root_writes_samples_only(tmp_path, command, output_tree, output_format):
    """`treetime` and `treetime clock` crashed in the root-to-tip plot and wrote an extra leaf (issue #959)."""
    from Bio import Phylo

    tree_file, dates_file = _write_single_child_root_inputs(tmp_path, SINGLE_CHILD_ROOT_DATES)
    outdir = tmp_path / 'out'

    _run_cli([*command, '--tree', tree_file, '--dates', dates_file, '--sequence-length', 1000, '--outdir', outdir])

    tree = Phylo.read(outdir / output_tree, output_format)
    assert sorted(n.name for n in tree.get_terminals()) == sorted(SINGLE_CHILD_ROOT_DATES)
    assert (outdir / 'root_to_tip_regression.pdf').is_file()


def test_cli_keep_root_keeps_single_child_root(tmp_path):
    from Bio import Phylo

    tree_file, dates_file = _write_single_child_root_inputs(tmp_path, SINGLE_CHILD_ROOT_DATES)
    outdir = tmp_path / 'out'

    _run_cli(['--tree', tree_file, '--dates', dates_file, '--sequence-length', 1000, '--keep-root', '--outdir', outdir])

    tree = Phylo.read(outdir / 'timetree.nexus', 'nexus')
    assert len(tree.root.clades) == 1
    assert sorted(n.name for n in tree.get_terminals()) == sorted(SINGLE_CHILD_ROOT_DATES)
    rows = [line.split('\t') for line in (outdir / 'dates.tsv').read_text().splitlines()[1:]]
    numeric_dates = {row[0]: float(row[2]) for row in rows}
    assert numeric_dates['ROOT'] < numeric_dates['MRCA']


def test_bad_branch_flags_are_python_bools():
    """Undated leaves G and H make their parent GH a bad branch."""
    dates = {name: date for name, date in SINGLE_CHILD_ROOT_DATES.items() if name not in ('G', 'H')}
    tt = _single_child_root_tree_time(dates)

    flags = {n.name: n.bad_branch for n in tt.tree.find_clades()}

    assert all(type(flag) is bool for flag in flags.values())
    assert {name for name, flag in flags.items() if flag} == {'G', 'H', 'GH'}


def test_tip_value_treats_falsy_bad_branch_as_good():
    import numpy as np

    tt = _single_child_root_tree_time(SINGLE_CHILD_ROOT_DATES)
    leaf = next(n for n in tt.tree.get_terminals() if n.name == 'A')
    leaf.bad_branch = np.False_

    assert tt.setup_TreeRegression().tip_value(leaf) == SINGLE_CHILD_ROOT_DATES['A']


def test_reroot_rejects_new_undated_leaf():
    """Without the single-child root removal, rerooting leaves ROOT as an undated leaf."""
    from unittest.mock import patch
    from treetime import TreeTime, TreeTimeError

    tt = _single_child_root_tree_time(SINGLE_CHILD_ROOT_DATES)

    with patch.object(TreeTime, '_remove_undated_single_child_root'):
        with pytest.raises(TreeTimeError, match='into leaves: ROOT$'):
            tt.reroot('least-squares')
