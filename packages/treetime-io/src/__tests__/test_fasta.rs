#[cfg(test)]
mod tests {
  use crate::fasta::{FastaRecord, fasta_read, fasta_write_record};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime_primitives::{AlphabetLike, AsciiChar, Seq};
  use treetime_utils::assert_error;

  #[test]
  fn test_fasta_reads_multiline_records_in_upper_case_with_descriptions() {
    let fasta = indoc! {"
      >a first sample
      acgt
      NN-a
      >b
      TTGC
    "};

    let actual = fasta_read(fasta.as_bytes(), &helpers::Nuc).unwrap();

    let expected = vec![
      helpers::record("a", Some("first sample"), "ACGTNN-A"),
      helpers::record("b", None, "TTGC"),
    ];
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_fasta_rejects_a_character_outside_the_alphabet() {
    let fasta = indoc! {"
      >a
      ACGT
      >b desc
      ACXT
    "};

    assert_error!(
      fasta_read(fasta.as_bytes(), &helpers::Nuc),
      r#"When processing sequence #2: ">b desc": FASTA input is incorrect: character "X" is not in the alphabet. Expected characters: '-', 'A', 'C', 'G', 'N', 'T'"#
    );
  }

  #[test]
  fn test_fasta_rejects_a_non_ascii_character() {
    let fasta = ">a\nAC\u{e9}T\n";

    assert_error!(
      fasta_read(fasta.as_bytes(), &helpers::Nuc),
      r#"When processing sequence #1: ">a": AsciiChar: 'é' is not ASCII"#
    );
  }

  #[test]
  fn test_fasta_write_record_writes_headers_with_descriptions() -> Result<(), Report> {
    let mut buf = Vec::new();
    fasta_write_record(&mut buf, "a", Some("first sample"), &Seq::try_from_str("ACGT")?)?;
    fasta_write_record(&mut buf, "b", None, &Seq::try_from_str("TT")?)?;

    assert_eq!(">a first sample\nACGT\n>b\nTT\n", String::from_utf8(buf)?);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) struct Nuc;

    impl AlphabetLike for Nuc {
      fn contains(&self, c: AsciiChar) -> bool {
        b"-ACGNT".contains(&u8::from(c))
      }

      fn chars(&self) -> impl Iterator<Item = AsciiChar> {
        b"-ACGNT".iter().map(|&byte| AsciiChar::from_byte_unchecked(byte))
      }
    }

    pub(super) fn record(name: &str, desc: Option<&str>, seq: &str) -> FastaRecord {
      FastaRecord {
        seq_name: name.to_owned(),
        desc: desc.map(str::to_owned),
        seq: Seq::try_from_str(seq).unwrap(),
      }
    }
  }
}
