"""
Test validation and error messages for SequenceData
"""
import pytest
from io import StringIO
from Bio import SeqIO
from treetime.sequence_data import SequenceData
from treetime import MissingDataError
import tempfile
import os


def test_inconsistent_sequence_lengths():
    """Test that mismatched sequence lengths produce a helpful error message"""

    # Create a temporary FASTA file with mismatched lengths
    fasta_content = """>seq1
ATGCATGC
>seq2
ATGCATGCATGC
>seq3
ATGCATGC
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        temp_file = f.name

    try:
        # This should raise a MissingDataError with a helpful message
        with pytest.raises(MissingDataError) as excinfo:
            sd = SequenceData(aln=temp_file)

        error_msg = str(excinfo.value)
        # Check that the error message contains key information
        assert "inconsistent lengths" in error_msg.lower()
        assert "8 positions" in error_msg or "12 positions" in error_msg
        assert "seq2" in error_msg  # The sequence with different length should be named

    finally:
        os.unlink(temp_file)


def test_duplicate_sequence_ids():
    """Test that duplicate sequence IDs produce a helpful error message"""

    # Create a temporary FASTA file with duplicate IDs
    fasta_content = """>seq1
ATGCATGC
>seq2
ATGCATGC
>seq1
GCTAGCTA
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        temp_file = f.name

    try:
        # This should raise a MissingDataError with a helpful message about duplicates
        with pytest.raises(MissingDataError) as excinfo:
            sd = SequenceData(aln=temp_file)

        error_msg = str(excinfo.value)
        # Check that the error message mentions duplicates
        assert "duplicate" in error_msg.lower()
        assert "seq1" in error_msg  # The duplicated sequence should be named

    finally:
        os.unlink(temp_file)


def test_valid_alignment_loads():
    """Test that a valid alignment loads without errors"""

    # Create a temporary FASTA file with valid sequences
    fasta_content = """>seq1
ATGCATGC
>seq2
ATGCATGC
>seq3
GCTAGCTA
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        temp_file = f.name

    try:
        # This should load successfully
        sd = SequenceData(aln=temp_file)
        assert sd is not None
        assert len(sd.aln) == 3

    finally:
        os.unlink(temp_file)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
