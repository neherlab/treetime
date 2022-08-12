# Test the arg functions on multiple trees
def test_reading_in_TK_output():
    import treetime.arg
    tree_nwk_files = ['test/arg/TreeKnit/tree_a_resolved.nwk','test/arg/TreeKnit/tree_b_resolved.nwk', 'test/arg/TreeKnit/tree_c_resolved.nwk']
    tree_names = ["tree_a", "tree_b", "tree_c"]
    assert treetime.arg.get_tree_names(tree_nwk_files) == tree_names
    MCC_dict = treetime.arg.get_MCC_dict('test/arg/TreeKnit/MCCs.json')

    # read trees, assert trees are named consistently
    tree_names = ["tree_a", "tree_b", "tree_c"]
    trees_in_dict = set().union(*MCC_dict.keys())
    assert all([k in trees_in_dict for k in tree_names])== True 

test_reading_in_TK_output()
def test_assign_mccs():
    from Bio import Phylo
    from treetime import TreeTime
    from treetime import arg
    from treetime.utils import parse_dates

    tree = Phylo.read('test/arg/TreeKnit/tree_a_resolved.nwk', 'newick')
    tt = TreeTime(dates=parse_dates("test/arg/TreeKnit/metadata.csv"), tree=tree,
            aln="test/arg/TreeKnit/aln_a.fasta", gtr='JC69', alphabet='nuc', verbose=True,
            fill_overhangs=True, keep_node_order=True,
            compress=False)

    # make a lookup for the MCCs and assign to tree
    MCC_dict = arg.get_MCC_dict('test/arg/TreeKnit/MCCs.json')
    MCC_locs = [frozenset(["tree_a", "tree_b"]), frozenset(["tree_a", "tree_c"])]
    MCCs = [MCC_dict[loc] for loc in MCC_locs]
    leaf_to_MCC = arg.get_mcc_map(MCCs)
    assert leaf_to_MCC == {'3_0': [0, 0], '10_0': [1, 1], '4_0': [1, 1], '5_0': [2, 2], '8_0': [2, 2], '1_0': [3, 3], '2_0': [3, 3], '6_0': [3, 3], '7_0': [3, 3], '9_0': [3, 3]}
    arg.assign_all_mccs(tt.tree, len(MCCs), leaf_to_MCC, tt.one_mutation)
    for node in tt.tree.find_clades():
        assert node.mcc[0] == node.mcc[1]
        if set([c.name for c in node.clades]) == set(["4_0", "10_0"]):
            node.mcc == [1,1]
        if set([c.name for c in node.clades]) == set(["8_0", "5_0"]):
            node.mcc == [2,2]

test_assign_mccs()

def test_parse_args():
    from treetime import arg
    import numpy as np
    
    tree_nwk_files = ['test/arg/TreeKnit/tree_a_resolved.nwk','test/arg/TreeKnit/tree_b_resolved.nwk', 'test/arg/TreeKnit/tree_c_resolved.nwk']
    aln_files = ['test/arg/TreeKnit/aln_a.fasta','test/arg/TreeKnit/aln_b.fasta', 'test/arg/TreeKnit/aln_c.fasta']
    MCC_file = 'test/arg/TreeKnit/MCCs.json'

    dict_ = arg.parse_arg(tree_nwk_files, aln_files, MCC_file, fill_overhangs=True)
    assert sum(dict_["masks_dict"][frozenset(["tree_a"])]) == sum(dict_["masks_dict"][frozenset(["tree_b"])]) == sum(dict_["masks_dict"][frozenset(["tree_c"])]) ==1000
    assert all(dict_["masks_dict"][frozenset(["tree_a"])] == np.concatenate((np.ones(1000), np.zeros(2000))))
    assert all(dict_["masks_dict"][frozenset(["tree_a", "tree_b"])] == np.concatenate((np.ones(2000), np.zeros(1000))))
    assert all(dict_["masks_dict"][frozenset(["tree_a", "tree_b", "tree_c"])] == np.ones(3000))

test_parse_args()

def test_setup_arg():
    from treetime import arg
    from treetime.utils import parse_dates
    import numpy as np

    tree_nwk_files = ['test/arg/TreeKnit/tree_a_resolved.nwk','test/arg/TreeKnit/tree_b_resolved.nwk', 'test/arg/TreeKnit/tree_c_resolved.nwk']
    aln_files = ['test/arg/TreeKnit/aln_a.fasta','test/arg/TreeKnit/aln_b.fasta', 'test/arg/TreeKnit/aln_c.fasta']
    MCC_file = 'test/arg/TreeKnit/MCCs.json'
    dates = parse_dates("test/arg/TreeKnit/metadata.csv")

    dict_ = arg.parse_arg(tree_nwk_files, aln_files, MCC_file, fill_overhangs=True)

    ##check if arg is set up correctly on tree_b
    masked_tree_b = arg.setup_arg(dict_["trees_dict"], dict_["alignment"], dates, dict_["MCCs_dict"], dict_["masks_dict"], "tree_b", gtr='JC69',
            verbose=0, fill_overhangs=True, reroot=False, fixed_clock_rate=0.001, alphabet='nuc')

    node_dict = {}
    for node in masked_tree_b.tree.find_clades():
        node_dict[node.name] = node
    for node in masked_tree_b.tree.find_clades():
        if node.name == "3_0" or set([c.name for c in node.clades]) == set(["4_0", "10_0"]) or set([c.name for c in node.clades]) == set(["3_0", "internal_9"]) or set([c.name for c in node.clades]) == set(["8_0", "5_0"]):
            assert all(node.mask == np.concatenate((np.zeros(1000), np.ones(2000))))
        elif node.name in set([c.name for c in masked_tree_b.tree.root.clades]):
            assert all(node.mask == np.concatenate((np.ones(2000), np.zeros(1000))))
        elif not node.up:
            assert all(node.mask == np.concatenate((np.zeros(1000), np.ones(1000), np.zeros(1000))))
        else:
            assert all(node.mask == np.ones(3000))

test_setup_arg()
