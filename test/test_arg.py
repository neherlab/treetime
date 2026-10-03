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
