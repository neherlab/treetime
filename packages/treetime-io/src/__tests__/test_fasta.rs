#[cfg(test)]
mod tests {
  use crate::fasta::{FastaRecord, read_many_fasta_str};
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

    let actual = read_many_fasta_str(fasta, &helpers::Nuc).unwrap();

    let expected = vec![
      helpers::record("a", Some("first sample"), "ACGTNN-A", 0),
      helpers::record("b", None, "TTGC", 1),
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
      read_many_fasta_str(fasta, &helpers::Nuc),
      r#"When processing sequence #2: ">b desc": FASTA input is incorrect: character "X" is not in the alphabet. Expected characters: '-', 'A', 'C', 'G', 'N', 'T'"#
    );
  }

  #[test]
  fn test_fasta_rejects_a_non_ascii_character() {
    let fasta = ">a\nAC\u{e9}T\n";

    assert_error!(
      read_many_fasta_str(fasta, &helpers::Nuc),
      r#"When processing sequence #1: ">a": AsciiChar: 'é' is not ASCII"#
    );
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

    pub(super) fn record(name: &str, desc: Option<&str>, seq: &str, index: usize) -> FastaRecord {
      FastaRecord {
        seq_name: name.to_owned(),
        desc: desc.map(str::to_owned),
        seq: Seq::try_from_str(seq).unwrap(),
        index,
      }
    }
  }
}
