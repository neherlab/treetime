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
    tt = TreeTime(dates=parse_dates("test/arg/metadata.csv"), tree=tree,
            aln="test/arg/aln_a.fasta", gtr='JC69', alphabet='nuc', verbose=True,
            fill_overhangs=True, keep_node_order=True,
            compress=False)

    # make a lookup for the MCCs and assign to tree
    MCC_dict = arg.get_MCC_dict('test/arg/TreeKnit/MCCs.json')
    MCC_locs = [frozenset(["tree_a", "tree_b"]), frozenset(["tree_a", "tree_c"])]
    MCCs = [MCC_dict[loc] for loc in MCC_locs]
    leaf_to_MCC = arg.get_mcc_map(MCCs)
    assert leaf_to_MCC == {'3': [0, 0], '6': [1, 1], '8': [2, 2], '4': [3, 3], '7': [3, 3], '1': [4, 3], '2': [4, 3], '5': [4, 3], '9': [4, 3], '10': [4, 3]}
    arg.assign_all_mccs(tt.tree, len(MCCs), leaf_to_MCC, tt.one_mutation)
    assert all([tt.tree.root.mcc == c.mcc for c in tt.tree.root.clades])
    assert sorted([c.mcc for c in tt.tree.root.clades[0]]) == sorted([[0, 0], [3,3], [4, 3], [4,3]])
    assert sorted([c.mcc for c in tt.tree.root.clades[1]]) == sorted([[1, 1], [4, 3]])
    assert sorted([c.mcc for c in tt.tree.root.clades[2]]) == sorted([[2, 2], [4, 3], [4,3]]) 

test_assign_mccs()

def test_parse_args():
    from treetime import arg
    import numpy as np
    
    tree_nwk_files = ['test/arg/TreeKnit/tree_a_resolved.nwk','test/arg/TreeKnit/tree_b_resolved.nwk', 'test/arg/TreeKnit/tree_c_resolved.nwk']
    aln_files = ['test/arg/aln_a.fasta','test/arg/aln_b.fasta', 'test/arg/aln_c.fasta']
    MCC_file = 'test/arg/TreeKnit/MCCs.json'

    dict_ = arg.parse_arg(tree_nwk_files, aln_files, MCC_file, fill_overhangs=True)
    assert sum(dict_["masks"][frozenset(["tree_a"])]) == sum(dict_["masks"][frozenset(["tree_b"])]) == sum(dict_["masks"][frozenset(["tree_c"])]) ==1000
    assert all(dict_["masks"][frozenset(["tree_a"])] == np.concatenate((np.ones(1000), np.zeros(2000))))
    assert all(dict_["masks"][frozenset(["tree_a", "tree_b"])] == np.concatenate((np.ones(2000), np.zeros(1000))))
    assert all(dict_["masks"][frozenset(["tree_a", "tree_b", "tree_c"])] == np.ones(3000))

test_parse_args()