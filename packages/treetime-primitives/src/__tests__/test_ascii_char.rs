#[cfg(test)]
mod tests {
  use crate::AsciiChar;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::lowest(  0)]
  #[case::highest( 127)]
  #[trace]
  fn test_ascii_char_try_new_accepts_ascii_range_bounds(#[case] value: u8) -> Result<(), Report> {
    let actual = AsciiChar::try_new(value)?;
    assert_eq!(value, actual.inner());
    Ok(())
  }

  #[test]
  fn test_ascii_char_try_new_rejects_first_non_ascii_byte() {
    let result = AsciiChar::try_new(128);
    assert_error!(result, "AsciiChar: value 128 is not ASCII (>= 128)");
  }

  #[test]
  fn test_ascii_char_try_from_char_keeps_ascii_byte() -> Result<(), Report> {
    let actual = AsciiChar::try_from_char('~')?;
    assert_eq!(b'~', actual.inner());
    Ok(())
  }

  #[test]
  fn test_ascii_char_try_from_char_rejects_non_ascii_char() {
    let result = AsciiChar::try_from_char('ñ');
    assert_error!(result, "AsciiChar: 'ñ' is not ASCII");
  }
}
