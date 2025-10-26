"""
Test tree loading with improved error messages
"""
import pytest
import tempfile
import os
from treetime import MissingDataError
from treetime.treeanc import read_tree


def test_invalid_newick_tree():
    """Test that invalid newick tree produces a helpful error message"""

    # Create a temporary file with invalid newick
    tree_content = """(A:0.1,B:0.2,C:0.3"""  # Missing closing paren

    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        f.write(tree_content)
        temp_file = f.name

    try:
        # This should raise a MissingDataError with details about what was tried
        with pytest.raises(MissingDataError) as excinfo:
            tree = read_tree(temp_file)

        error_msg = str(excinfo.value)
        # Check that the error message shows what was tried
        assert "failed to read tree" in error_msg.lower()
        assert "attempted formats" in error_msg.lower()
        # Should show both newick and nexus attempts
        assert "newick" in error_msg.lower()
        assert "nexus" in error_msg.lower()

    finally:
        os.unlink(temp_file)


def test_valid_newick_tree():
    """Test that a valid newick tree loads successfully"""

    # Create a temporary file with valid newick
    tree_content = """(A:0.1,B:0.2,C:0.3);"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        f.write(tree_content)
        temp_file = f.name

    try:
        # This should load successfully
        tree = read_tree(temp_file)
        assert tree is not None
        assert tree.count_terminals() == 3

    finally:
        os.unlink(temp_file)


def test_valid_nexus_tree():
    """Test that a valid nexus tree loads successfully"""

    # Create a temporary file with valid nexus
    nexus_content = """#NEXUS
begin trees;
    tree tree1 = (A:0.1,B:0.2,C:0.3);
end;
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.nex', delete=False) as f:
        f.write(nexus_content)
        temp_file = f.name

    try:
        # This should load successfully
        tree = read_tree(temp_file)
        assert tree is not None
        assert tree.count_terminals() == 3

    finally:
        os.unlink(temp_file)


def test_tree_too_few_terminals():
    """Test that trees with too few terminals produce a deprecation warning"""

    # Create a tree with only 2 terminals (minimum is 3)
    tree_content = """(A:0.1,B:0.2);"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
        f.write(tree_content)
        temp_file = f.name

    try:
        # This should emit a DeprecationWarning but load successfully
        with pytest.warns(DeprecationWarning) as warning_list:
            tree = read_tree(temp_file)

        # Check that the warning message mentions terminals and future error
        assert len(warning_list) == 1
        warning_msg = str(warning_list[0].message)
        assert "terminals" in warning_msg.lower()
        assert "0.12" in warning_msg  # Should mention the version

        # Check that the tree loaded successfully
        assert tree is not None
        assert tree.count_terminals() == 2

    finally:
        os.unlink(temp_file)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
