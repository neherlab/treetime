use eyre::Report;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
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

    let position_str = std::str::from_utf8(&bytes[1..bytes.len() - 1])
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
  use treetime_utils::assert_error;

  #[test]
  fn test_mutation_parse_a123t() -> Result<(), Report> {
    let result = Mutation::parse("A123T")?;
    let expected = Mutation::new(b'A', 123, b'T');
    assert_eq!(expected, result);
    Ok(())
  }

  #[test]
  fn test_mutation_parse_c1g() -> Result<(), Report> {
    let result = Mutation::parse("C1G")?;
    let expected = Mutation::new(b'C', 1, b'G');
    assert_eq!(expected, result);
    Ok(())
  }

  #[test]
  fn test_mutation_parse_t99999a() -> Result<(), Report> {
    let result = Mutation::parse("T99999A")?;
    let expected = Mutation::new(b'T', 99999, b'A');
    assert_eq!(expected, result);
    Ok(())
  }

  #[test]
  fn test_mutation_parse_lowercase() -> Result<(), Report> {
    let result = Mutation::parse("a1b")?;
    let expected = Mutation::new(b'a', 1, b'b');
    assert_eq!(expected, result);
    Ok(())
  }

  #[test]
  fn test_mutation_parse_with_whitespace() -> Result<(), Report> {
    let result = Mutation::parse(" G42C ")?;
    let expected = Mutation::new(b'G', 42, b'C');
    assert_eq!(expected, result);
    Ok(())
  }

  #[test]
  fn test_mutation_parse_empty_fails() {
    assert_error!(Mutation::parse(""), "Mutation string too short: ''");
  }

  #[test]
  fn test_mutation_parse_single_char_fails() {
    assert_error!(Mutation::parse("A"), "Mutation string too short: 'A'");
  }

  #[test]
  fn test_mutation_parse_two_chars_fails() {
    assert_error!(Mutation::parse("AT"), "Mutation string too short: 'AT'");
  }

  #[test]
  fn test_mutation_parse_only_digits_fails() {
    assert_error!(Mutation::parse("123"), "Invalid reference character in mutation: '123'");
  }

  #[test]
  fn test_mutation_parse_digit_first_fails() {
    assert_error!(Mutation::parse("1A2"), "Invalid reference character in mutation: '1A2'");
  }

  #[test]
  fn test_mutation_parse_negative_position_fails() {
    assert_error!(
      Mutation::parse("A-1T"),
      "Invalid position in mutation 'A-1T': invalid digit found in string"
    );
  }

  #[test]
  fn test_mutation_format_a123t() {
    let mutation = Mutation::new(b'A', 123, b'T');
    assert_eq!("A123T", mutation.format());
  }

  #[test]
  fn test_mutation_format_c1g() {
    let mutation = Mutation::new(b'C', 1, b'G');
    assert_eq!("C1G", mutation.format());
  }

  #[test]
  fn test_mutation_format_t99999a() {
    let mutation = Mutation::new(b'T', 99999, b'A');
    assert_eq!("T99999A", mutation.format());
  }

  #[test]
  fn test_mutation_roundtrip() -> Result<(), Report> {
    let inputs = ["A123T", "C1G", "T99999A", "G42C"];
    for input in inputs {
      let mutation = Mutation::parse(input)?;
      let formatted = mutation.format();
      assert_eq!(input, formatted);
    }
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
