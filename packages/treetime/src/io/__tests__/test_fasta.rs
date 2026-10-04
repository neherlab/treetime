#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::o;
  use eyre::Report;
  use helpers::{record, seq};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::sync::LazyLock;
  use treetime_io::fasta::{FastaRecord, fasta_read};
  use treetime_utils::assert_error;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);
  static AA_ALPHABET: LazyLock<Alphabet> = LazyLock::new(|| Alphabet::new(AlphabetName::Aa).unwrap());

  #[test]
  fn test_fasta_reader_fail_on_non_fasta() {
    let data =
        b"This is not a valid FASTA string.\nIt is not empty, and not entirely whitespace\nbut does not contain 'greater than' character.\n";
    assert_error!(
      fasta_read(data.as_slice(), &*NUC_ALPHABET),
      "FASTA input is incorrectly formatted: expected at least one FASTA record starting with character '>', but none found"
    );
  }

  #[test]
  fn test_fasta_reader_fail_on_unknown_char() {
    let data = b">seq%1\nACGT%ACGT\n";
    assert_error!(
      fasta_read(data.as_slice(), &*NUC_ALPHABET),
      r#"When processing sequence #1: ">seq%1": FASTA input is incorrect: character "%" is not in the alphabet. Expected characters: '-', 'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'V', 'W', 'Y'"#
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(                       b"".as_slice(),                              vec![])]
  #[case::whitespace_only(             b"\n \n \n\n".as_slice(),                    vec![])]
  #[case::single_record(               b">seq1\nATCG\n".as_slice(),                 vec![record("seq1", "ATCG")])]
  #[case::leading_newline(             b"\n>seq1\nATCG\n".as_slice(),               vec![record("seq1", "ATCG")])]
  #[case::multiple_leading_newlines(   b"\n\n\n>seq1\nATCG\n".as_slice(),           vec![record("seq1", "ATCG")])]
  #[case::no_trailing_newline(         b">seq1\nATCG".as_slice(),                   vec![record("seq1", "ATCG")])]
  #[case::trailing_empty_line(         b">seq1\nATCG\n\n".as_slice(),               vec![record("seq1", "ATCG")])]
  #[case::multiple_records(            b">seq1\nATCG\n>seq2\nGCTA\n".as_slice(),    vec![record("seq1", "ATCG"), record("seq2", "GCTA")])]
  #[case::empty_lines_between_records( b"\n>seq1\n\nATCG\n\n\n>seq2\nGCTA\n\n".as_slice(), vec![record("seq1", "ATCG"), record("seq2", "GCTA")])]
  #[case::leading_newlines_last_unterminated(b"\n\n>a\nACGCTCGATC\n\n>b\nCCGCGC".as_slice(), vec![record("a", "ACGCTCGATC"), record("b", "CCGCGC")])]
  #[case::last_record_without_sequence(b">a\nACGCTCGATC\n>b\nCCGCGC\n>c".as_slice(), vec![record("a", "ACGCTCGATC"), record("b", "CCGCGC"), record("c", "")])]
  #[case::middle_record_without_sequence(b">a\nACGCTCGATC\n>b\n>c\nCCGCGC".as_slice(), vec![record("a", "ACGCTCGATC"), record("b", ""), record("c", "CCGCGC")])]
  #[case::first_record_empty(          b">\n>C\nACGT\n>D\nACGA\n".as_slice(),      vec![record("", ""), record("C", "ACGT"), record("D", "ACGA")])]
  #[trace]
  fn test_fasta_reader_records(#[case] data: &[u8], #[case] expected: Vec<FastaRecord>) -> Result<(), Report> {
    let actual = fasta_read(data, &*NUC_ALPHABET)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_name_desc() -> Result<(), Report> {
    let actual = fasta_read(
      indoc! {r#"
        >Identifier Description
        ACGT
        >Identifier Description with spaces
        ACGT


      "#}
      .as_bytes(),
      &*NUC_ALPHABET,
    )?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("Identifier"),
        desc: Some(o!("Description")),
        seq: seq("ACGT"),
      },
      FastaRecord {
        seq_name: o!("Identifier"),
        desc: Some(o!("Description with spaces")),
        seq: seq("ACGT"),
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  #[cfg_attr(
    dylint_lib = "treetime_lints",
    expect(
      typographic_characters,
      reason = "sequence names with emoji exercise Unicode handling"
    )
  )]
  fn test_fasta_reader_dedent_nuc() -> Result<(), Report> {
    let actual = fasta_read(
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
      "#}
      .as_bytes(),
      &*NUC_ALPHABET,
    )?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("FluBuster-001"),
        desc: None,
        seq: seq("ACAGCCATGTATTG--"),
      },
      FastaRecord {
        seq_name: o!("CommonCold-AB"),
        desc: None,
        seq: seq("ACATCCCTGTA-TG--"),
      },
      FastaRecord {
        seq_name: o!("Ecoli/Joke/2024|XD"),
        desc: None,
        seq: seq("ACATCGCCNNA--GAC"),
      },
      FastaRecord {
        seq_name: o!("Sniffles-B"),
        desc: None,
        seq: seq("GCATCCCTGTA-NG--"),
      },
      FastaRecord {
        seq_name: o!("StrawberryYogurtCulture|🍓"),
        desc: None,
        seq: seq("CCGGCCATGTATTG--"),
      },
      FastaRecord {
        seq_name: o!(""),
        desc: Some(o!("SneezeC-19")),
        seq: seq("CCGGCGATGTRTTG--"),
      },
      FastaRecord {
        seq_name: o!("MisindentedVirus|D-skew"),
        desc: None,
        seq: seq("TCGGCCGTGTRTTG--"),
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  #[cfg_attr(
    dylint_lib = "treetime_lints",
    expect(
      typographic_characters,
      reason = "sequence names with emoji exercise Unicode handling"
    )
  )]
  fn test_fasta_reader_dedent_aa() -> Result<(), Report> {
    let actual = fasta_read(
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
      "#}
      .as_bytes(),
      &*AA_ALPHABET,
    )?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("Prot/000|β-Napkinase"),
        desc: None,
        seq: seq("MXDXXXTQ-B--"),
      },
      FastaRecord {
        seq_name: o!("Enzyme/2024|LaughzymeFactor"),
        desc: None,
        seq: seq("AX*XB-TQVWR*"),
      },
      FastaRecord {
        seq_name: o!("😊-Gigglecatalyst"),
        desc: None,
        seq: seq("MKXTQWX-B**"),
      },
      FastaRecord {
        seq_name: o!("CellFunSignal"),
        desc: None,
        seq: seq("MQXQXXBQRW**"),
      },
      FastaRecord {
        seq_name: o!("Pathway/042|Doodlease"),
        desc: None,
        seq: seq("MXQ-*XTQWBQR"),
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fasta_reader_multiline_and_skewed_indentation() -> Result<(), Report> {
    let actual = fasta_read(
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
      "#}
      .as_bytes(),
      &*NUC_ALPHABET,
    )?;

    let expected = vec![
      FastaRecord {
        seq_name: o!("MixedCaseSeq"),
        desc: None,
        seq: seq("ACAGCCATGTATTG--"),
      },
      FastaRecord {
        seq_name: o!("LowercaseSeq"),
        desc: None,
        seq: seq("ACAGCCATGTATTG--"),
      },
      FastaRecord {
        seq_name: o!("UppercaseSeq"),
        desc: None,
        seq: seq("ACAGCCATGTATTG--"),
      },
      FastaRecord {
        seq_name: o!("MultilineSeq"),
        desc: None,
        seq: seq("ACAGCCATGTATTG--"),
      },
      FastaRecord {
        seq_name: o!("SkewedIndentSeq"),
        desc: None,
        seq: seq("ACAGCCATGTATTGATTG--"),
      },
    ];

    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::o;
    use treetime_io::fasta::FastaRecord;
    use treetime_primitives::Seq;

    pub(super) fn seq(s: &str) -> Seq {
      Seq::try_from_str(s).unwrap()
    }

    pub(super) fn record(name: &str, sequence: &str) -> FastaRecord {
      FastaRecord {
        seq_name: o!(name),
        desc: None,
        seq: seq(sequence),
      }
    }
  }
}
