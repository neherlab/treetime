use crate::seq_char::AsciiChar;
use eyre::Report;
use treetime_utils::error::make_error;

/// Represents genetic sequence (ASCII characters only)
#[must_use]
#[derive(Clone, PartialOrd, Ord, Default)]
pub struct Seq {
  data: Vec<AsciiChar>,
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

  pub fn try_from_elem<T: Into<u8>>(elem: T, n: usize) -> Result<Self, Report> {
    let ch = AsciiChar::try_new(elem.into())?;
    Ok(Self { data: vec![ch; n] })
  }

  /// Create a sequence from an ASCII string.
  #[allow(unsafe_code)]
  pub fn try_from_str(s: &str) -> Result<Self, Report> {
    if !s.is_ascii() {
      return make_error!("Seq: input contains non-ASCII characters");
    }
    // SAFETY: We just validated that s.is_ascii() above
    Ok(unsafe { Self::from_str_unchecked(s) })
  }

  /// Create a sequence from a pre-validated ASCII string.
  ///
  /// # Safety
  ///
  /// The caller must ensure that `s` contains only ASCII characters (bytes 0-127).
  /// Passing non-ASCII input violates the type invariant and causes undefined behavior
  /// when calling [`as_str`](Self::as_str).
  #[allow(unsafe_code)]
  pub unsafe fn from_str_unchecked(s: &str) -> Self {
    debug_assert!(
      s.is_ascii(),
      "Seq::from_str_unchecked: input contains non-ASCII characters"
    );
    Self {
      data: s
        .as_bytes()
        .iter()
        .copied()
        .map(AsciiChar::from_byte_unchecked)
        .collect(),
    }
  }

  /// Create a sequence from a vector of bytes.
  pub fn try_from_vec(vec: Vec<u8>) -> Result<Self, Report> {
    let data = vec.into_iter().map(AsciiChar::try_new).collect::<Result<Vec<_>, _>>()?;
    Ok(Self { data })
  }

  /// Create a sequence from a byte slice.
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

  pub fn truncate(&mut self, len: usize) {
    self.data.truncate(len);
  }

  pub fn push(&mut self, byte: AsciiChar) {
    self.data.push(byte);
  }

  pub fn pop(&mut self) -> Option<AsciiChar> {
    self.data.pop()
  }

  pub fn capacity(&self) -> usize {
    self.data.capacity()
  }

  pub fn reserve(&mut self, additional: usize) {
    self.data.reserve(additional);
  }

  pub fn reserve_exact(&mut self, additional: usize) {
    self.data.reserve_exact(additional);
  }

  #[allow(unsafe_code)]
  pub fn as_str(&self) -> &str {
    // SAFETY: `self.data.as_ptr()` is guaranteed to be valid for reads and properly aligned
    // because `data` is a `Vec<AsciiChar>`. The `AsciiChar` type ensures that each element is a valid
    // single-byte ASCII character, making the conversion to a `u8` pointer valid.
    let byte_slice = unsafe { std::slice::from_raw_parts(self.data.as_ptr().cast::<u8>(), self.data.len()) };

    // SAFETY: `from_utf8_unchecked` is safe here because `byte_slice` is guaranteed to contain only
    // valid UTF-8 data. This is ensured by the invariant that `AsciiChar` can only hold valid ASCII characters,
    // which are a subset of UTF-8.
    unsafe { std::str::from_utf8_unchecked(byte_slice) }
  }

  pub fn as_slice(&self) -> &[AsciiChar] {
    &self.data
  }

  pub fn as_mut_slice(&mut self) -> &mut [AsciiChar] {
    &mut self.data
  }

  pub fn try_insert(&mut self, index: usize, byte: u8) -> Result<(), Report> {
    let ch = AsciiChar::try_new(byte)?;
    self.data.insert(index, ch);
    Ok(())
  }

  pub fn remove(&mut self, index: usize) -> u8 {
    self.data.remove(index).into()
  }

  pub fn try_append(&mut self, other: &mut Vec<u8>) -> Result<(), Report> {
    for byte in other.drain(..) {
      self.data.push(AsciiChar::try_new(byte)?);
    }
    Ok(())
  }

  pub fn try_push_str(&mut self, s: &str) -> Result<(), Report> {
    if !s.is_ascii() {
      return make_error!("Seq: input contains non-ASCII characters");
    }
    self
      .data
      .extend(s.as_bytes().iter().copied().map(AsciiChar::from_byte_unchecked));
    Ok(())
  }

  pub fn contains_str(&self, s: &str) -> bool {
    self.as_str().contains(s)
  }

  pub fn starts_with(&self, prefix: &str) -> bool {
    self.as_str().starts_with(prefix)
  }

  pub fn ends_with(&self, suffix: &str) -> bool {
    self.as_str().ends_with(suffix)
  }

  pub fn find(&self, substring: &str) -> Option<usize> {
    self.as_str().find(substring)
  }

  pub fn rfind(&self, substring: &str) -> Option<usize> {
    self.as_str().rfind(substring)
  }

  pub fn try_replace(&mut self, from: &str, to: &str) -> Result<&mut Self, Report> {
    let replaced = self.as_str().replace(from, to);
    self.data = replaced
      .into_bytes()
      .into_iter()
      .map(AsciiChar::try_new)
      .collect::<Result<Vec<_>, _>>()?;
    Ok(self)
  }

  pub fn split_off(&mut self, at: usize) -> Seq {
    Seq {
      data: self.data.split_off(at),
    }
  }

  // pub fn splice<I>(
  //   &mut self,
  //   range: core::ops::Range<usize>,
  //   replace_with: I,
  // ) -> std::vec::Splice<'_, <I as IntoIterator>::IntoIter>
  // where
  //   I: IntoIterator<Item = u8>,
  // {
  //   self.data.splice(range, replace_with.into_iter().map(AsciiChar::from))
  // }
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

impl core::fmt::Display for Seq {
  fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
    f.write_str(self.as_str())
  }
}

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

impl IntoIterator for Seq {
  type Item = AsciiChar;
  type IntoIter = std::vec::IntoIter<AsciiChar>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.into_iter()
  }
}

#[allow(unsafe_code)]
impl std::io::Read for Seq {
  fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
    let len = std::cmp::min(buf.len(), self.len());
    // SAFETY:
    // 1. `self.data` is guaranteed to hold only ASCII characters because `Seq` enforces this invariant via `AsciiChar`.
    // 2. The length `len` is calculated as the minimum of `buf.len()` and `self.len()`, ensuring no out-of-bounds access for either slice.
    // 3. `std::ptr::copy_nonoverlapping` is safe to use here because:
    //    a. Both `self.data` and `buf` are valid, properly aligned, and non-overlapping.
    //    b. The memory regions are guaranteed to be accessible for `len` bytes.
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

#[allow(unsafe_code)]
impl<'de> serde::Deserialize<'de> for Seq {
  fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
  where
    D: serde::Deserializer<'de>,
  {
    let s = String::deserialize(deserializer)?;
    if !s.is_ascii() {
      return Err(serde::de::Error::custom("Seq: input contains non-ASCII characters"));
    }
    // SAFETY: We just validated that s.is_ascii() above
    Ok(unsafe { Seq::from_str_unchecked(&s) })
  }
}

/// Create a `Seq` from `AsciiChar` values.
///
/// - `seq!()` creates an empty sequence
/// - `seq![ch; n]` creates a sequence with `n` copies of `AsciiChar` `ch`
/// - `seq![ch1, ch2, ...]` creates a sequence from `AsciiChar` values
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
