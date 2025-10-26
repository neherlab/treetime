"""Test validation and error messages for SequenceData"""
import pytest
from treetime.sequence_data import SequenceData
from treetime import MissingDataError
import tempfile
import os


def test_alignment_error_messages():
    """Test that alignment errors preserve BioPython messages and show all format attempts"""
    # Mismatched sequence lengths
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">seq1\nATGCATGC\n>seq2\nATGCATGCATGC\n")
        temp_file = f.name

    try:
        with pytest.raises(MissingDataError) as exc:
            SequenceData(aln=temp_file)
        assert "attempted formats" in str(exc.value).lower()
        assert "same length" in str(exc.value).lower()
    finally:
        os.unlink(temp_file)


def test_duplicate_ids_warning():
    """Test that duplicate IDs emit deprecation warning"""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(">seq1\nATGC\n>seq1\nGCTA\n")
        temp_file = f.name

    try:
        with pytest.warns(DeprecationWarning, match="(?si)duplicate.*seq1.*0\\.12"):
            sd = SequenceData(aln=temp_file)
        assert len(sd.aln) == 1  # BioPython keeps only last occurrence
    finally:
        os.unlink(temp_file)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
