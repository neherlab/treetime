#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum NhxType {
  Decimal,
  Integer,
  Duplication,
  Color,
  Text,
}

impl NhxType {
  pub(crate) const fn of(key: &str) -> Self {
    match key.as_bytes() {
      b"B" => Self::Decimal,
      b"T" | b"W" => Self::Integer,
      b"D" => Self::Duplication,
      b"C" => Self::Color,
      _ => Self::Text,
    }
  }

  pub(crate) const fn description(self) -> &'static str {
    match self {
      Self::Decimal => "a decimal number",
      Self::Integer => "an integer",
      Self::Duplication => "one of T, F, Y, N or ?",
      Self::Color => "a color written as red.green.blue",
      Self::Text => "text",
    }
  }
}

pub(crate) const DUPLICATION_VALUES: [&str; 5] = ["T", "F", "Y", "N", "?"];
