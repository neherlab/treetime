use crate::seq_char::AsciiChar;
use eyre::Report;
use treetime_utils::error::make_error;

impl<'a> IntoIterator for &'a Seq {
  type Item = &'a AsciiChar;
  type IntoIter = core::slice::Iter<'a, AsciiChar>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.iter()
  }
}

impl<'a> IntoIterator for &'a mut Seq {
  type Item = &'a mut AsciiChar;
  type IntoIter = core::slice::IterMut<'a, AsciiChar>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.iter_mut()
  }
}

#[must_use]
#[derive(Clone, PartialOrd, Ord, Default)]
pub struct Seq {
  data: Vec<AsciiChar>,
}

impl PartialEq<Seq> for str {
  fn eq(&self, other: &Seq) -> bool {
    other == self
  }
}

impl PartialEq<Seq> for String {
  fn eq(&self, other: &Seq) -> bool {
    other == self.as_str()
  }
}

impl Seq {
  pub fn new() -> Self {
    Self { data: Vec::new() }
  }

  pub fn with_capacity(capacity: usize) -> Self {
    Self {
      data: Vec::with_capacity(capacity),
    }
  }

  pub fn try_from_str(s: &str) -> Result<Self, Report> {
    if !s.is_ascii() {
      return make_error!("Seq: input contains non-ASCII characters");
    }
    Ok(Self::from_ascii_str(s))
  }

  fn from_ascii_str(s: &str) -> Self {
    debug_assert!(s.is_ascii(), "Seq::from_ascii_str: input contains non-ASCII characters");
    Self {
      data: s
        .as_bytes()
        .iter()
        .copied()
        .map(AsciiChar::from_byte_unchecked)
        .collect(),
    }
  }

  pub fn try_from_slice(slice: &[u8]) -> Result<Self, Report> {
    let data = slice
      .iter()
      .copied()
      .map(AsciiChar::try_new)
      .collect::<Result<Vec<_>, _>>()?;
    Ok(Self { data })
  }

  pub fn len(&self) -> usize {
    self.data.len()
  }

  pub fn is_empty(&self) -> bool {
    self.data.is_empty()
  }

  pub fn clear(&mut self) {
    self.data.clear();
  }

  pub fn push(&mut self, byte: AsciiChar) {
    self.data.push(byte);
  }

  pub fn reserve(&mut self, additional: usize) {
    self.data.reserve(additional);
  }

  #[allow(
    unsafe_code,
    clippy::undocumented_unsafe_blocks,
    reason = "AsciiChar is transparent over a validated ASCII byte"
  )]
  pub fn as_str(&self) -> &str {
    let byte_slice = unsafe { std::slice::from_raw_parts(self.data.as_ptr().cast::<u8>(), self.data.len()) };

    unsafe { std::str::from_utf8_unchecked(byte_slice) }
  }

  pub fn as_slice(&self) -> &[AsciiChar] {
    &self.data
  }

  pub fn as_mut_slice(&mut self) -> &mut [AsciiChar] {
    &mut self.data
  }
}

impl PartialEq for Seq {
  fn eq(&self, other: &Self) -> bool {
    self.data == other.data
  }
}

impl PartialEq<Vec<u8>> for Seq {
  fn eq(&self, other: &Vec<u8>) -> bool {
    if self.data.len() != other.len() {
      return false;
    }
    self.data.iter().zip(other).all(|(c, &b)| c.inner() == b)
  }
}

impl PartialEq<str> for Seq {
  fn eq(&self, other: &str) -> bool {
    if self.data.len() != other.len() {
      return false;
    }
    self.data.iter().zip(other.as_bytes()).all(|(c, &b)| c.inner() == b)
  }
}

impl PartialEq<String> for Seq {
  fn eq(&self, other: &String) -> bool {
    self == other.as_str()
  }
}

impl PartialEq<&str> for Seq {
  fn eq(&self, other: &&str) -> bool {
    self == *other
  }
}

impl Eq for Seq {}

impl core::ops::Deref for Seq {
  type Target = [AsciiChar];
  fn deref(&self) -> &Self::Target {
    &self.data
  }
}

impl core::ops::DerefMut for Seq {
  fn deref_mut(&mut self) -> &mut Self::Target {
    &mut self.data
  }
}

impl From<&[AsciiChar]> for Seq {
  fn from(slice: &[AsciiChar]) -> Self {
    slice.iter().copied().collect()
  }
}

impl Extend<AsciiChar> for Seq {
  fn extend<I: IntoIterator<Item = AsciiChar>>(&mut self, iter: I) {
    self.data.extend(iter);
  }
}

impl FromIterator<AsciiChar> for Seq {
  fn from_iter<I: IntoIterator<Item = AsciiChar>>(iter: I) -> Self {
    Self {
      data: Vec::from_iter(iter),
    }
  }
}

impl AsRef<[AsciiChar]> for Seq {
  fn as_ref(&self) -> &[AsciiChar] {
    &self.data
  }
}

impl AsRef<[u8]> for Seq {
  fn as_ref(&self) -> &[u8] {
    self.as_str().as_bytes()
  }
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(handwritten_fmt_impl, reason = "a sequence renders as its characters")
)]
impl core::fmt::Display for Seq {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    f.write_str(self.as_str())
  }
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    handwritten_fmt_impl,
    reason = "Debug shows the sequence text, not its byte vector, for readable test diffs"
  )
)]
impl core::fmt::Debug for Seq {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    core::fmt::Display::fmt(self, f)
  }
}

impl core::ops::Index<usize> for Seq {
  type Output = AsciiChar;

  fn index(&self, index: usize) -> &Self::Output {
    &self.data[index]
  }
}

impl core::ops::IndexMut<usize> for Seq {
  fn index_mut(&mut self, index: usize) -> &mut Self::Output {
    &mut self.data[index]
  }
}

impl core::ops::Index<core::ops::Range<usize>> for Seq {
  type Output = [AsciiChar];

  fn index(&self, index: core::ops::Range<usize>) -> &Self::Output {
    &self.data[index]
  }
}

impl core::ops::IndexMut<core::ops::Range<usize>> for Seq {
  fn index_mut(&mut self, index: core::ops::Range<usize>) -> &mut Self::Output {
    &mut self.data[index]
  }
}

impl std::ops::Add for Seq {
  type Output = Self;

  fn add(mut self, other: Self) -> Self::Output {
    self.data.extend(other.data);
    self
  }
}

impl std::ops::Mul<usize> for Seq {
  type Output = Self;

  fn mul(mut self, rhs: usize) -> Self::Output {
    let original = self.data.clone();
    for _ in 1..rhs {
      self.data.extend(&original);
    }
    self
  }
}

impl IntoIterator for Seq {
  type Item = AsciiChar;
  type IntoIter = std::vec::IntoIter<AsciiChar>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.into_iter()
  }
}

#[allow(
  unsafe_code,
  clippy::undocumented_unsafe_blocks,
  reason = "AsciiChar is transparent over a validated ASCII byte"
)]
impl std::io::Read for Seq {
  fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
    let len = std::cmp::min(buf.len(), self.len());
    unsafe {
      std::ptr::copy_nonoverlapping(self.data.as_ptr().cast::<u8>(), buf.as_mut_ptr(), len);
    }
    self.data.drain(..len);
    Ok(len)
  }
}

impl std::io::Write for Seq {
  fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
    for &byte in buf {
      let ch =
        AsciiChar::try_new(byte).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))?;
      self.data.push(ch);
    }
    Ok(buf.len())
  }

  fn flush(&mut self) -> std::io::Result<()> {
    Ok(())
  }
}

impl serde::Serialize for Seq {
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: serde::Serializer,
  {
    serializer
      .serialize_str(self.as_str())
      .map_err(serde::ser::Error::custom)
  }
}

impl<'de> serde::Deserialize<'de> for Seq {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: serde::Deserializer<'de>,
  {
    let s = String::deserialize(deserializer)?;
    if !s.is_ascii() {
      return Err(serde::de::Error::custom("Seq: input contains non-ASCII characters"));
    }
    Ok(Seq::from_ascii_str(&s))
  }
}

#[macro_export]
macro_rules! seq {
  () => (
      $crate::seq::Seq::new()
  );
  ($elem:expr; $n:expr) => (
      $crate::seq::Seq::from_iter(::std::iter::repeat_n($elem, $n))
  );
  ($($char:expr),* $(,)?) => {
    {
      $crate::seq::Seq::from_iter([$($char),*].into_iter())
    }
  };
}
