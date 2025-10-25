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


def test_phylip_inconsistent_lengths():
    """Test that phylip files with mismatched lengths produce helpful error"""

    # Create a phylip file with mismatched lengths
    phylip_content = """3 8
seq1      ATGCATGC
seq2      ATGCATGCATGC
seq3      ATGCATGC
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.phy', delete=False) as f:
        f.write(phylip_content)
        temp_file = f.name

    try:
        with pytest.raises(MissingDataError) as excinfo:
            sd = SequenceData(aln=temp_file)

        error_msg = str(excinfo.value)
        # Should either get our detailed message or the original BioPython error
        assert "length" in error_msg.lower()

    finally:
        os.unlink(temp_file)


def test_nexus_inconsistent_lengths():
    """Test that nexus files with mismatched lengths produce helpful error"""

    # Create a nexus file with mismatched lengths
    nexus_content = """#NEXUS
begin data;
dimensions ntax=3 nchar=8;
format datatype=dna;
matrix
seq1 ATGCATGC
seq2 ATGCATGCATGC
seq3 ATGCATGC
;
end;
"""

    with tempfile.NamedTemporaryFile(mode='w', suffix='.nex', delete=False) as f:
        f.write(nexus_content)
        temp_file = f.name

    try:
        with pytest.raises(MissingDataError) as excinfo:
            sd = SequenceData(aln=temp_file)

        error_msg = str(excinfo.value)
        # Should catch the nexus error about data length
        assert "loading alignment failed" in error_msg.lower()
        assert ("nchar" in error_msg.lower() or "length" in error_msg.lower())

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
