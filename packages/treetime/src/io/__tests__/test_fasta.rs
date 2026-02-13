#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::o;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use pretty_assertions::assert_eq;
  use std::io::Cursor;
  use treetime_io::fasta::*;
  use treetime_utils::error::report_to_string;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
    static ref AA_ALPHABET: Alphabet = Alphabet::new(AlphabetName::Aa, false).unwrap();
  }

  #[test]
  fn test_fasta_reader_fail_on_non_fasta() {
    let data =
        b"This is not a valid FASTA string.\nIt is not empty, and not entirely whitespace\nbut does not contain 'greater than' character.\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);
    let mut record = FastaRecord::new();
    assert_eq!(
      reader.read(&mut record).unwrap_err().to_string(),
      "FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found"
    );
  }

  #[test]
  fn test_fasta_reader_fail_on_unknown_char() {
    let data = b">seq%1\nACGT%ACGT\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);
    let mut record = FastaRecord::new();
    let actual = report_to_string(&reader.read(&mut record).unwrap_err());
    let expected = r#"When processing sequence #1: ">seq%1": FASTA input is incorrect: character "%" is not in the alphabet. Expected characters: '-', 'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'V', 'W', 'Y'"#;
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_fasta_reader_read_empty() {
    let data = b"";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert!(record.is_empty());
  }

  #[test]
  fn test_fasta_reader_read_whitespace_only() {
    let data = b"\n \n \n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert!(record.is_empty());
  }

  #[test]
  fn test_fasta_reader_read_single_record() {
    let data = b">seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_with_leading_newline() {
    let data = b"\n>seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_with_multiple_leading_newlines() {
    let data = b"\n\n\n>seq1\nATCG\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_single_record_without_trailing_newline() {
    let data = b">seq1\nATCG";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_read_multiple_records() {
    let data = b">seq1\nATCG\n>seq2\nGCTA\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record1 = FastaRecord::new();
    reader.read(&mut record1).unwrap();

    let mut record2 = FastaRecord::new();
    reader.read(&mut record2).unwrap();

    assert_eq!(record1.seq_name, "seq1");
    assert_eq!(record1.seq, "ATCG");
    assert_eq!(record1.index, 0);

    assert_eq!(record2.seq_name, "seq2");
    assert_eq!(record2.seq, "GCTA");
    assert_eq!(record2.index, 1);
  }

  #[test]
  fn test_fasta_reader_read_empty_lines_between_records() {
    let data = b"\n>seq1\n\nATCG\n\n\n>seq2\nGCTA\n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record1 = FastaRecord::new();
    reader.read(&mut record1).unwrap();

    let mut record2 = FastaRecord::new();
    reader.read(&mut record2).unwrap();

    assert_eq!(record1.seq_name, "seq1");
    assert_eq!(record1.seq, "ATCG");
    assert_eq!(record1.index, 0);

    assert_eq!(record2.seq_name, "seq2");
    assert_eq!(record2.seq, "GCTA");
    assert_eq!(record2.index, 1);
  }

  #[test]
  fn test_fasta_reader_read_with_trailing_newline() {
    let data = b">seq1\nATCG\n\n";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(record.seq_name, "seq1");
    assert_eq!(record.seq, "ATCG");
    assert_eq!(record.index, 0);
  }

  #[test]
  fn test_fasta_reader_example_1() {
    let data = b"\n\n>a\nACGCTCGATC\n\n>b\nCCGCGC";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: "ACGCTCGATC".into(),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: "CCGCGC".into(),
        index: 1,
      }
    );
  }

  #[test]
  fn test_fasta_reader_example_2() {
    let data = b">a\nACGCTCGATC\n>b\nCCGCGC\n>c";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: "ACGCTCGATC".into(),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: "CCGCGC".into(),
        index: 1,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("c"),
        desc: None,
        seq: "".into(),
        index: 2,
      }
    );
  }

  #[test]
  fn test_fasta_reader_example_3() {
    let data = b">a\nACGCTCGATC\n>b\n>c\nCCGCGC";
    let mut reader = FastaReader::new(Box::new(Cursor::new(data)), &*NUC_ALPHABET);

    let mut record = FastaRecord::new();
    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("a"),
        desc: None,
        seq: "ACGCTCGATC".into(),
        index: 0,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("b"),
        desc: None,
        seq: "".into(),
        index: 1,
      }
    );

    reader.read(&mut record).unwrap();

    assert_eq!(
      record,
      FastaRecord {
        seq_name: o!("c"),
        desc: None,
        seq: "CCGCGC".into(),
        index: 2,
      }
    );
  }

  #[test]
  fn test_fasta_reader_name_desc() -> Result<(), Report> {
    let actual = read_many_fasta_str(
      indoc! {r#"
        >Identifier Description
        ACGT
        >Identifier Description with spaces
        ACGT


      "#},
      &*NUC_ALPHABET,
    )?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("Identifier"),
        desc: Some(o!("Description")),
        seq: "ACGT".into(),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("Identifier"),
        desc: Some(o!("Description with spaces")),
        seq: "ACGT".into(),
        index: 1,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_dedent_nuc() -> Result<(), Report> {
    let actual = read_many_fasta_str(
      indoc! {r#"
        >FluBuster-001
        ACAGCCATGTATTG--
        >CommonCold-AB
        ACATCCCTGTA-TG--
        >Ecoli/Joke/2024|XD
        ACATCGCCNNA--GAC

        >Sniffles-B
        GCATCCCTGTA-NG--
        >StrawberryYogurtCulture|🍓
        CCGGCCATGTATTG--
        > SneezeC-19
        CCGGCGATGTRTTG--
          >MisindentedVirus|D-skew
          TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("FluBuster-001"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("CommonCold-AB"),
        desc: None,
        seq: "ACATCCCTGTA-TG--".into(),
        index: 1,
      },
      FastaRecord {
        seq_name: o!("Ecoli/Joke/2024|XD"),
        desc: None,
        seq: "ACATCGCCNNA--GAC".into(),
        index: 2,
      },
      FastaRecord {
        seq_name: o!("Sniffles-B"),
        desc: None,
        seq: "GCATCCCTGTA-NG--".into(),
        index: 3,
      },
      FastaRecord {
        seq_name: o!("StrawberryYogurtCulture|🍓"),
        desc: None,
        seq: "CCGGCCATGTATTG--".into(),
        index: 4,
      },
      FastaRecord {
        seq_name: o!(""),
        desc: Some(o!("SneezeC-19")),
        seq: "CCGGCGATGTRTTG--".into(),
        index: 5,
      },
      FastaRecord {
        seq_name: o!("MisindentedVirus|D-skew"),
        desc: None,
        seq: "TCGGCCGTGTRTTG--".into(),
        index: 6,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_dedent_aa() -> Result<(), Report> {
    let actual = read_many_fasta_str(
      indoc! {r#"
        >Prot/000|β-Napkinase
        MXDXXXTQ-B--
        >Enzyme/2024|LaughzymeFactor
        AX*XB-TQVWR*

        >😊-Gigglecatalyst
        MKXTQWX-B**
        >CellFunSignal
        MQXQXXBQRW**
        >Pathway/042|Doodlease
        MXQ-*XTQWBQR
      "#},
      &*AA_ALPHABET,
    )?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("Prot/000|β-Napkinase"),
        desc: None,
        seq: "MXDXXXTQ-B--".into(),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("Enzyme/2024|LaughzymeFactor"),
        desc: None,
        seq: "AX*XB-TQVWR*".into(),
        index: 1,
      },
      FastaRecord {
        seq_name: o!("😊-Gigglecatalyst"),
        desc: None,
        seq: "MKXTQWX-B**".into(),
        index: 2,
      },
      FastaRecord {
        seq_name: o!("CellFunSignal"),
        desc: None,
        seq: "MQXQXXBQRW**".into(),
        index: 3,
      },
      FastaRecord {
        seq_name: o!("Pathway/042|Doodlease"),
        desc: None,
        seq: "MXQ-*XTQWBQR".into(),
        index: 4,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_multiline_and_skewed_indentation() -> Result<(), Report> {
    let actual = read_many_fasta_str(
      indoc! {r#"
        >MixedCaseSeq
        aCaGcCAtGtAtTG--
        >LowercaseSeq
        acagccatgtattg--
        >UppercaseSeq
        ACAGCCATGTATTG--
        >MultilineSeq
        ACAGCC
        ATGT
        ATTG--
        >SkewedIndentSeq
          ACAGCC
        ATGTATTG
         ATTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("MixedCaseSeq"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 0,
      },
      FastaRecord {
        seq_name: o!("LowercaseSeq"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 1,
      },
      FastaRecord {
        seq_name: o!("UppercaseSeq"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 2,
      },
      FastaRecord {
        seq_name: o!("MultilineSeq"),
        desc: None,
        seq: "ACAGCCATGTATTG--".into(),
        index: 3,
      },
      FastaRecord {
        seq_name: o!("SkewedIndentSeq"),
        desc: None,
        seq: "ACAGCCATGTATTGATTG--".into(),
        index: 4,
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }
}
