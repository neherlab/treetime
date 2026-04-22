use eyre::Report;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::str::from_utf8;
use treetime::make_error;

pub type MutationList = Vec<Mutation>;
pub type PartitionedMutations = BTreeMap<String, MutationList>;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct Mutation {
  pub position: usize,
  pub reference: u8,
  pub alternative: u8,
}

impl Mutation {
  pub fn new(reference: u8, position: usize, alternative: u8) -> Self {
    Self {
      position,
      reference,
      alternative,
    }
  }

  pub fn parse(s: &str) -> Result<Self, Report> {
    let s = s.trim();
    if s.len() < 3 {
      return make_error!("Mutation string too short: '{s}'");
    }

    let bytes = s.as_bytes();
    let reference = bytes[0];
    let alternative = bytes[bytes.len() - 1];

    if !reference.is_ascii_alphabetic() {
      return make_error!("Invalid reference character in mutation: '{s}'");
    }
    if !alternative.is_ascii_alphabetic() {
      return make_error!("Invalid alternative character in mutation: '{s}'");
    }

    let position_str = from_utf8(&bytes[1..bytes.len() - 1])
      .map_err(|e| treetime::make_report!("Invalid UTF-8 in mutation '{s}': {e}"))?;
    let position = position_str
      .parse::<usize>()
      .map_err(|e| treetime::make_report!("Invalid position in mutation '{s}': {e}"))?;

    Ok(Self {
      position,
      reference,
      alternative,
    })
  }

  pub fn format(&self) -> String {
    format!(
      "{}{}{}",
      self.reference as char, self.position, self.alternative as char
    )
  }
}

pub fn parse_mutation_list(mutations: &[String]) -> Result<MutationList, Report> {
  mutations.iter().map(|s| Mutation::parse(s)).try_collect()
}

pub fn format_mutation_list(mutations: &[Mutation]) -> Vec<String> {
  mutations.iter().map(Mutation::format).collect_vec()
}

#[cfg(test)]
mod tests {
  use super::*;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rstest]
  #[case::standard("A123T", b'A', 123, b'T')]
  #[case::single_digit_pos("C1G", b'C', 1, b'G')]
  #[case::large_position("T99999A", b'T', 99999, b'A')]
  #[case::lowercase("a1b", b'a', 1, b'b')]
  #[case::with_whitespace(" G42C ", b'G', 42, b'C')]
  fn test_mutation_parse(
    #[case] input: &str,
    #[case] reference: u8,
    #[case] position: usize,
    #[case] alternative: u8,
  ) -> Result<(), Report> {
    let expected = Mutation::new(reference, position, alternative);
    let result = Mutation::parse(input)?;
    assert_eq!(expected, result);
    Ok(())
  }

  #[rstest]
  #[case::empty("", "Mutation string too short: ''")]
  #[case::single_char("A", "Mutation string too short: 'A'")]
  #[case::two_chars("AT", "Mutation string too short: 'AT'")]
  #[case::no_letters("123", "Invalid reference character in mutation: '123'")]
  #[case::digit_first("1A2", "Invalid reference character in mutation: '1A2'")]
  #[case::negative_position("A-1T", "Invalid position in mutation 'A-1T': invalid digit found in string")]
  fn test_mutation_parse_error(#[case] input: &str, #[case] expected_error: &str) {
    assert_error!(Mutation::parse(input), expected_error);
  }

  #[rstest]
  #[case::standard(b'A', 123, b'T', "A123T")]
  #[case::single_digit_pos(b'C', 1, b'G', "C1G")]
  #[case::large_position(b'T', 99999, b'A', "T99999A")]
  fn test_mutation_format(
    #[case] reference: u8,
    #[case] position: usize,
    #[case] alternative: u8,
    #[case] expected: &str,
  ) {
    let mutation = Mutation::new(reference, position, alternative);
    assert_eq!(expected, mutation.format());
  }

  #[rstest]
  #[case::standard("A123T")]
  #[case::single_digit_pos("C1G")]
  #[case::large_position("T99999A")]
  #[case::two_digit_pos("G42C")]
  fn test_mutation_roundtrip(#[case] input: &str) -> Result<(), Report> {
    let mutation = Mutation::parse(input)?;
    let formatted = mutation.format();
    assert_eq!(input, formatted);
    Ok(())
  }

  #[test]
  fn test_parse_mutation_list() -> Result<(), Report> {
    let input = vec!["A1T".to_owned(), "C2G".to_owned(), "G3A".to_owned()];
    let expected = vec![
      Mutation::new(b'A', 1, b'T'),
      Mutation::new(b'C', 2, b'G'),
      Mutation::new(b'G', 3, b'A'),
    ];
    let result = parse_mutation_list(&input)?;
    assert_eq!(expected, result);
    Ok(())
  }

  #[test]
  fn test_format_mutation_list() {
    let input = vec![
      Mutation::new(b'A', 1, b'T'),
      Mutation::new(b'C', 2, b'G'),
      Mutation::new(b'G', 3, b'A'),
    ];
    let expected = vec!["A1T".to_owned(), "C2G".to_owned(), "G3A".to_owned()];
    let result = format_mutation_list(&input);
    assert_eq!(expected, result);
  }
}
