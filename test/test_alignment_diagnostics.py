import pytest
from treetime import MissingDataError
from treetime.sequence_data import SequenceData


def test_mismatched_fasta_lengths_reports_offending_sequences(tmp_path):
    fasta = tmp_path / 'bad_lengths.fasta'
    seq_171 = 'M' * 171
    fasta.write_text(
        f'>seq1\n{seq_171}\n'
        f'>seq2\n{seq_171}\n'
        f'>seq3\n{seq_171}\n'
        f'>AP019183\n{"A" * 722}\n'
        f'>E01011\n{"C" * 751}\n'
        f'>FJ428830\n{"G" * 751}\n',
        encoding='utf-8',
    )
    with pytest.raises(MissingDataError) as exc_info:
        SequenceData(str(fasta))

    msg = str(exc_info.value)
    assert 'SequenceData: loading alignment failed' in msg
    assert 'Expected: 171 characters (3 of 6 sequences)' in msg
    assert 'AP019183: 722 characters' in msg
    assert 'E01011: 751 characters' in msg
    assert 'FJ428830: 751 characters' in msg
    assert 'All sequences in an alignment must have the same length.' in msg


def test_empty_alignment_file_reports_empty_diagnostic(tmp_path):
    empty_fasta = tmp_path / 'empty.fasta'
    empty_fasta.write_text('', encoding='utf-8')
    with pytest.raises(MissingDataError) as exc_info:
        SequenceData(str(empty_fasta))

    msg = str(exc_info.value)
    assert 'SequenceData: loading alignment failed' in msg
    assert 'is empty' in msg


def test_missing_alignment_file_reports_missing_diagnostic(tmp_path):
    missing = tmp_path / 'does_not_exist.fasta'
    with pytest.raises(MissingDataError) as exc_info:
        SequenceData(str(missing))

    msg = str(exc_info.value)
    assert 'SequenceData: loading alignment failed' in msg
    assert 'does not exist' in msg


def test_duplicate_ids_with_length_mismatch_noted(tmp_path):
    fasta = tmp_path / 'dup_and_mismatch.fasta'
    fasta.write_text(
        '>MT648657\nAAAA\n>MT648657\nAAAA\n>OUTLIER\nAAAAAA\n',
        encoding='utf-8',
    )
    with pytest.raises(MissingDataError) as exc_info:
        SequenceData(str(fasta))

    msg = str(exc_info.value)
    assert 'OUTLIER: 6 characters' in msg
    assert 'Duplicate sequence IDs' in msg
    assert 'MT648657' in msg
