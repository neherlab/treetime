#[cfg(test)]
mod tests {
  use crate::Seq;

  #[test]
  #[should_panic(expected = "Seq::from_str: input contains non-ASCII characters")]
  fn test_from_str_non_ascii_panics() {
    let _unused = Seq::from_str("hello\u{0080}world");
  }

  #[test]
  #[should_panic(expected = "AsciiChar::from(u8): value 128 >= 128")]
  fn test_from_vec_non_ascii_panics() {
    let _unused = Seq::from_vec(vec![65, 66, 128, 67]);
  }

  #[test]
  #[should_panic(expected = "AsciiChar::from(u8): value 200 >= 128")]
  fn test_from_slice_non_ascii_panics() {
    let _unused = Seq::from_slice(&[65, 66, 200, 67]);
  }
}
