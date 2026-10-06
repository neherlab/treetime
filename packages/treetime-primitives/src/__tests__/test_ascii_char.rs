#[cfg(test)]
mod tests {
  use crate::AsciiChar;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::json_read_str;

  #[test]
  fn test_try_new_valid() -> Result<(), Report> {
    let actual = AsciiChar::try_new(127)?;
    let expected = AsciiChar::from_byte_unchecked(b'~' + 1);
    assert_eq!(actual.inner(), expected.inner());
    Ok(())
  }

  #[test]
  fn test_try_new_invalid_error() {
    let result = AsciiChar::try_new(128);
    assert_error!(result, "AsciiChar: value 128 is not ASCII (>= 128)");
  }

  #[test]
  fn test_try_new_zero_valid() -> Result<(), Report> {
    let actual = AsciiChar::try_new(0)?;
    assert_eq!(actual.inner(), 0);
    Ok(())
  }

  #[test]
  fn test_try_from_char_valid() -> Result<(), Report> {
    let actual = AsciiChar::try_from_char('~')?;
    assert_eq!(actual.inner(), b'~');
    Ok(())
  }

  #[test]
  fn test_try_from_char_invalid_error() {
    let result = AsciiChar::try_from_char('ñ');
    assert_error!(result, "AsciiChar: 'ñ' is not ASCII");
  }

  #[test]
  fn test_deserialize_valid() -> Result<(), Report> {
    let actual: AsciiChar = json_read_str("65")?;
    let expected = AsciiChar::from_byte_unchecked(b'A');
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_deserialize_invalid_fails() {
    assert_error!(
      json_read_str::<AsciiChar>("200"),
      "When parsing JSON: Unexpected: invalid value: AsciiChar: value 200 is not ASCII (>= 128) at line 1 column 1"
    );
  }
}
