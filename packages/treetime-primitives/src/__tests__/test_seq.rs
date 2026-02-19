#[cfg(test)]
mod tests {
  use crate::Seq;
  use treetime_utils::assert_error;

  #[test]
  fn test_try_from_str_non_ascii_error() {
    let result = Seq::try_from_str("hello\u{0080}world");
    assert_error!(result, "Seq: input contains non-ASCII characters");
  }

  #[test]
  fn test_try_from_vec_non_ascii_error() {
    let result = Seq::try_from_vec(vec![65, 66, 128, 67]);
    assert_error!(result, "AsciiChar: value 128 is not ASCII (>= 128)");
  }

  #[test]
  fn test_try_from_slice_non_ascii_error() {
    let result = Seq::try_from_slice(&[65, 66, 200, 67]);
    assert_error!(result, "AsciiChar: value 200 is not ASCII (>= 128)");
  }
}
