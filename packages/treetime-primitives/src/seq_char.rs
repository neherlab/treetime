use deser::adapters::{DeserializeAs, TryFromInto};
use deser::de::SinkHandle;
use deser::{Deserialize, Serialize, State};
use eyre::Report;
use std::fmt::Write as _;
use treetime_utils::error::make_error;

impl From<AsciiChar> for u8 {
  fn from(item: AsciiChar) -> Self {
    item.0
  }
}

impl From<AsciiChar> for u16 {
  fn from(item: AsciiChar) -> Self {
    u16::from(item.0)
  }
}

impl From<AsciiChar> for u32 {
  fn from(item: AsciiChar) -> Self {
    u32::from(item.0)
  }
}

impl From<AsciiChar> for u64 {
  fn from(item: AsciiChar) -> Self {
    u64::from(item.0)
  }
}

impl From<AsciiChar> for usize {
  fn from(item: AsciiChar) -> Self {
    usize::from(item.0)
  }
}

impl From<AsciiChar> for char {
  fn from(item: AsciiChar) -> Self {
    char::from(item.0)
  }
}

#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize)]
#[repr(transparent)]
pub struct AsciiChar(u8);

impl AsciiChar {
  pub fn try_new(value: u8) -> Result<Self, Report> {
    if value >= 128 {
      return make_error!("AsciiChar: value {value} is not ASCII (>= 128)");
    }
    Ok(Self(value))
  }

  #[allow(
    clippy::as_conversions,
    reason = "narrowing char to u8 is exact after the is_ascii guard"
  )]
  pub fn try_from_char(value: char) -> Result<Self, Report> {
    if !value.is_ascii() {
      return make_error!("AsciiChar: '{value}' is not ASCII");
    }
    Ok(Self(value as u8))
  }

  pub const fn from_byte_unchecked(value: u8) -> Self {
    debug_assert!(value < 128, "AsciiChar::from_byte_unchecked: value >= 128");
    Self(value)
  }

  pub const fn inner(&self) -> u8 {
    self.0
  }
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    handwritten_fmt_impl,
    reason = "a sequence character renders as the character, not its byte value"
  )
)]
impl core::fmt::Display for AsciiChar {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    f.write_char(char::from(self.0))
  }
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    handwritten_fmt_impl,
    reason = "Debug shows the character, not its byte value, for readable test diffs"
  )
)]
impl core::fmt::Debug for AsciiChar {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    core::fmt::Display::fmt(self, f)
  }
}

impl<'de> Deserialize<'de> for AsciiChar {
  fn deserialize_into<'out>(out: &'out mut Option<Self>, state: &mut State) -> SinkHandle<'out, 'de> {
    <TryFromInto<u8> as DeserializeAs<'de, Self>>::deserialize_into_as(out, state)
  }
}

impl TryFrom<u8> for AsciiChar {
  type Error = Report;

  fn try_from(value: u8) -> Result<Self, Report> {
    Self::try_new(value)
  }
}
