use std::fmt::Display;

#[macro_export]
macro_rules! o {
  ($x:expr $(,)?) => {
    ToOwned::to_owned($x)
  };
}

pub fn vec_to_string(v: Vec<char>) -> String {
  // Surprisingly, this is the fastest way, according to `benches/vec_char_to_string.rs`
  let bytes: Vec<u8> = v.into_iter().map(|c| c as u8).collect();
  String::from_utf8(bytes).unwrap()
}

pub fn quote(x: impl Display) -> String {
  format!("\"{x}\"")
}

pub fn quote_single(x: impl Display) -> String {
  format!("'{x}'")
}
