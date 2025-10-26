"""Test tree loading with improved error messages"""
import pytest
import tempfile
import os
from treetime import MissingDataError
from treetime.treeanc import read_tree


def test_tree_error_messages():
    """Test that tree errors show all format attempts with BioPython messages"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        f.write("(A:0.1,B:0.2,C:0.3")  # Invalid: missing closing paren
        temp_file = f.name

    try:
        with pytest.raises(MissingDataError) as exc:
            read_tree(temp_file)
        err = str(exc.value).lower()
        assert "attempted formats" in err
        assert "newick" in err and "nexus" in err
    finally:
        os.unlink(temp_file)


def test_tree_few_terminals_warning():
    """Test that trees with <3 terminals emit deprecation warning"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        f.write("(A:0.1,B:0.2);")
        temp_file = f.name

    try:
        with pytest.warns(DeprecationWarning, match="(?s)terminals.*0\\.12"):
            tree = read_tree(temp_file)
        assert tree.count_terminals() == 2
    finally:
        os.unlink(temp_file)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
