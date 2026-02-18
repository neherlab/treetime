#[cfg(test)]
mod tests {
  use crate::AsciiChar;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_new_valid() {
    let actual = AsciiChar::new(127);
    let expected = AsciiChar::new(b'~' + 1);
    assert_eq!(actual.inner(), expected.inner());
  }

  #[test]
  #[should_panic(expected = "AsciiChar::new: value >= 128")]
  fn test_new_invalid_panics() {
    let _ = AsciiChar::new(128);
  }

  #[test]
  fn test_from_u8_valid() {
    let actual = AsciiChar::from(0_u8);
    assert_eq!(actual.inner(), 0);
  }

  #[test]
  #[should_panic(expected = "AsciiChar::from(u8): value 128 >= 128")]
  fn test_from_u8_invalid_panics() {
    let _ = AsciiChar::from(128_u8);
  }

  #[test]
  fn test_from_char_valid() {
    let actual = AsciiChar::from('~');
    assert_eq!(actual.inner(), b'~');
  }

  #[test]
  #[should_panic(expected = "AsciiChar::from(char): 'ñ' is not ASCII")]
  fn test_from_char_invalid_panics() {
    let _ = AsciiChar::from('ñ');
  }

  #[test]
  fn test_deserialize_valid() {
    let actual: AsciiChar = serde_json::from_str("65").unwrap();
    let expected = AsciiChar::new(b'A');
    assert_eq!(actual, expected);
  }

  #[test]
  fn test_deserialize_invalid_fails() {
    let result: Result<AsciiChar, _> = serde_json::from_str("200");
    assert!(result.is_err());
    let err = result.unwrap_err().to_string();
    assert!(err.contains("200"), "error should mention the value: {err}");
    assert!(err.contains("128"), "error should mention the limit: {err}");
  }
}
