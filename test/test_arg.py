import json

import pytest

from treetime.arg import read_mccs, tree_label

MCCS = [['A'], ['B', 'C'], ['D', 'E', 'F']]


def write(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return str(path)


def test_read_mccs_text(tmp_path):
    path = write(tmp_path, 'MCCs.dat', 'A\nB,C\nD,E,F')
    assert read_mccs(path) == MCCS


def test_read_mccs_json_single_pair(tmp_path):
    content = json.dumps({'MCC_dict': {'1': {'trees': ['ha', 'na'], 'mccs': MCCS}}})
    path = write(tmp_path, 'MCCs.json', content)
    # a single pair is used whatever the tree labels
    assert read_mccs(path) == MCCS
    assert read_mccs(path, tree_labels=['tree1', 'tree2']) == MCCS


def test_read_mccs_json_selects_pair(tmp_path):
    other = [['A', 'B', 'C', 'D', 'E', 'F']]
    content = json.dumps(
        {
            'MCC_dict': {
                '1': {'trees': ['ha', 'na'], 'mccs': other},
                '2': {'trees': ['ha', 'pb2'], 'mccs': MCCS},
            }
        }
    )
    path = write(tmp_path, 'MCCs.json', content)
    assert read_mccs(path, tree_labels=['pb2', 'ha']) == MCCS
    with pytest.raises(ValueError, match='several tree pairs'):
        read_mccs(path, tree_labels=['ha', 'mp'])


def test_tree_label():
    assert tree_label('results/tree_ha.nwk') == 'tree_ha'
    # output trees of TreeKnit carry the label of their input tree
    assert tree_label('results/tree_ha_resolved.nwk') == 'tree_ha'
    assert tree_label('results/ARG/tree_na_liberal_resolved.nwk') == 'tree_na'
    assert tree_label('results/tree_pb2_imputed.nwk') == 'tree_pb2'


def phylo(newick):
    from io import StringIO

    from Bio import Phylo

    t = Phylo.read(StringIO(newick), 'newick')
    t.root.up = None
    for n in t.find_clades(order='preorder'):
        for c in n:
            c.up = n
    return t


def node(tree, *leaves):
    return tree.common_ancestor(*leaves)


def test_map_mccs():
    from treetime.arg import map_mccs

    t = phylo('((A,B),(C,(D,X)));')
    m = map_mccs(t, {'X': 0, 'A': 1, 'B': 1, 'C': 1, 'D': 1})
    assert m[node(t, 'D', 'X')] == 1  # parent's MCC continues into it
    assert m[node(t, 'A', 'B')] == 1
    assert m[t.root] == 1
    assert m[node(t, 'X')] == 0
    # leaves without an MCC (absent from the other segment) constrain nothing
    m = map_mccs(t, {'A': 1, 'B': 1, 'C': 1, 'D': 1})
    assert m[node(t, 'D', 'X')] == 1 and m[node(t, 'X')] is None


def test_masks_three_segments():
    import numpy as np

    from treetime.arg import annotate_mccs, assign_masks

    t = phylo('((A,B),(C,(D,X)));')
    seg = [np.array(x, dtype=float) for x in ([1, 0, 0], [0, 1, 0], [0, 0, 1])]
    # segment b shares everything with the focal tree; in segment c, X reassorted
    annotate_mccs(t, {'b': [['A', 'B', 'C', 'D', 'X']], 'c': [['X'], ['A', 'B', 'C', 'D']]})
    assign_masks(t, seg[0], {'b': seg[1], 'c': seg[2]})
    assert list(node(t, 'X').mask) == [1, 1, 0]
    assert list(node(t, 'D').mask) == [1, 1, 1]
    assert list(node(t, 'D', 'X').mask) == [1, 1, 1]
    assert list(t.root.mask) == [1, 0, 0]  # no branch above the root
    assert node(t, 'X').mccs == {'b': 0, 'c': 0}
    assert node(t, 'X').mcc is None  # only set with a single other segment


def test_mcc_comment():
    from treetime.CLI_io import mcc_comment

    t = phylo('(A,B);')
    n = node(t, 'A')
    n.mccs, n.mcc = {'na': 2}, 2
    assert mcc_comment(n) == ',mcc="2"'
    n.mccs, n.mcc = {'na': 2, 'pb2': None}, None
    assert mcc_comment(n) == ',mcc_na="2",mcc_pb2="None"'


def test_parse_arg_three_trees(tmp_path):
    from treetime.arg import parse_arg

    trees = []
    alns = []
    for k, (nwk, seqs) in enumerate(
        [
            ('((A,B),C);', {'A': 'AC', 'B': 'AG', 'C': 'TT'}),
            ('((A,C),B);', {'A': 'AAA', 'B': 'CCC', 'C': 'GGG'}),
            ('(A,(B,C));', {'A': 'T', 'B': 'T'}),  # C has no sequence for this segment
        ]
    ):
        trees.append(write(tmp_path, f'seg{k}.nwk', nwk))
        alns.append(write(tmp_path, f'seg{k}.fasta', ''.join(f'>{n}\n{s}\n' for n, s in seqs.items())))
    pairs = {
        '1': {'trees': ['seg0', 'seg1'], 'mccs': [['B'], ['A', 'C']]},
        '2': {'trees': ['seg0', 'seg2'], 'mccs': [['A'], ['B', 'C']]},
        '3': {'trees': ['seg1', 'seg2'], 'mccs': [['A', 'B', 'C']]},
    }
    mccs = write(tmp_path, 'MCCs.json', json.dumps({'MCC_dict': pairs}))
    p = parse_arg(trees, alns, mccs)
    assert p['labels'] == ['seg0', 'seg1', 'seg2']
    assert p['MCCs'] == {(0, 1): [['B'], ['A', 'C']], (0, 2): [['A'], ['B', 'C']], (1, 2): [['A', 'B', 'C']]}
    seqs = {s.id: str(s.seq) for s in p['alignment']}
    assert seqs == {'A': 'ACAAAT', 'B': 'AGCCCT', 'C': 'TTGGGN'}
    assert [list(m) for m in p['masks']] == [[1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 1, 0], [0, 0, 0, 0, 0, 1]]
