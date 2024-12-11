/// Represents genetic sequence (ASCII characters only)
#[must_use]
#[derive(Clone, PartialOrd, Ord, Default)]
pub struct Seq {
  data: Vec<u8>,
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

  pub fn from_elem<T: Into<u8>>(elem: T, n: usize) -> Self {
    Self {
      data: vec![elem.into(); n],
    }
  }

  pub fn from_string(s: String) -> Self {
    debug_assert!(s.is_ascii());
    Self { data: s.into_bytes() }
  }

  pub fn from_str(s: &str) -> Self {
    debug_assert!(s.is_ascii());
    Self {
      data: s.as_bytes().to_vec(),
    }
  }

  pub fn from_vec(vec: Vec<u8>) -> Self {
    Self { data: vec }
  }

  pub fn from_slice(slice: &[u8]) -> Self {
    Self { data: slice.to_vec() }
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

  pub fn push(&mut self, byte: u8) {
    self.data.push(byte);
  }

  pub fn pop(&mut self) -> Option<u8> {
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
    // SAFETY: safe if `data` consists of ASCII-characters
    unsafe { core::str::from_utf8_unchecked(&self.data) }
  }

  pub fn as_slice(&self) -> &[u8] {
    &self.data
  }

  pub fn as_mut_slice(&mut self) -> &mut [u8] {
    &mut self.data
  }

  pub fn insert(&mut self, index: usize, byte: u8) {
    self.data.insert(index, byte);
  }

  pub fn remove(&mut self, index: usize) -> u8 {
    self.data.remove(index)
  }

  pub fn append(&mut self, other: &mut Vec<u8>) {
    self.data.append(other);
  }

  pub fn push_str(&mut self, s: &str) {
    self.data.extend_from_slice(s.as_bytes());
  }

  pub fn contains(&self, s: &str) -> bool {
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

  pub fn replace(&mut self, from: &str, to: &str) -> &mut Self {
    let replaced = self.as_str().replace(from, to);
    self.data = replaced.into_bytes();
    self
  }

  pub fn split_off(&mut self, at: usize) -> Seq {
    Seq {
      data: self.data.split_off(at),
    }
  }

  pub fn splice<I>(
    &mut self,
    range: core::ops::Range<usize>,
    replace_with: I,
  ) -> std::vec::Splice<'_, <I as IntoIterator>::IntoIter>
  where
    I: IntoIterator<Item = u8>,
  {
    self.data.splice(range, replace_with)
  }
}

impl PartialEq for Seq {
  fn eq(&self, other: &Self) -> bool {
    self.data == other.data
  }
}

impl PartialEq<Vec<u8>> for Seq {
  fn eq(&self, other: &Vec<u8>) -> bool {
    self.data == *other
  }
}

impl PartialEq<str> for Seq {
  fn eq(&self, other: &str) -> bool {
    self.data == other.as_bytes()
  }
}

impl PartialEq<String> for Seq {
  fn eq(&self, other: &String) -> bool {
    self.data == other.as_bytes()
  }
}

impl PartialEq<&str> for Seq {
  fn eq(&self, other: &&str) -> bool {
    self.data == other.as_bytes()
  }
}

impl PartialEq<Seq> for str {
  fn eq(&self, other: &Seq) -> bool {
    self.as_bytes() == other.data
  }
}

impl PartialEq<Seq> for String {
  fn eq(&self, other: &Seq) -> bool {
    self.as_bytes() == other.data
  }
}

impl Eq for Seq {}

impl core::ops::Deref for Seq {
  type Target = [u8];
  fn deref(&self) -> &Self::Target {
    &self.data
  }
}

impl core::ops::DerefMut for Seq {
  fn deref_mut(&mut self) -> &mut Self::Target {
    &mut self.data
  }
}

impl From<&str> for Seq {
  fn from(s: &str) -> Self {
    Self::from_str(s)
  }
}

impl From<&[u8]> for Seq {
  fn from(slice: &[u8]) -> Self {
    Self { data: slice.to_vec() }
  }
}

impl Extend<u8> for Seq {
  fn extend<I: IntoIterator<Item = u8>>(&mut self, iter: I) {
    self.data.extend(iter);
  }
}

impl FromIterator<u8> for Seq {
  fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
    Self {
      data: Vec::from_iter(iter),
    }
  }
}

impl AsRef<[u8]> for Seq {
  fn as_ref(&self) -> &[u8] {
    &self.data
  }
}

impl AsRef<str> for Seq {
  fn as_ref(&self) -> &str {
    self.as_str()
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
  type Output = u8;

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
  type Output = [u8];

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
  type Item = &'a u8;
  type IntoIter = core::slice::Iter<'a, u8>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.iter()
  }
}

impl<'a> IntoIterator for &'a mut Seq {
  type Item = &'a mut u8;
  type IntoIter = core::slice::IterMut<'a, u8>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.iter_mut()
  }
}

impl IntoIterator for Seq {
  type Item = u8;
  type IntoIter = std::vec::IntoIter<u8>;

  fn into_iter(self) -> Self::IntoIter {
    self.data.into_iter()
  }
}

impl std::io::Write for Seq {
  fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
    self.data.extend_from_slice(buf);
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
    Ok(Seq::from_string(s))
  }
}

#[macro_export]
macro_rules! seq {
  () => (
      $crate::representation::seq::Seq::new()
  );
  ($elem:expr; $n:expr) => (
      $crate::representation::seq::Seq::from_elem($elem, $n)
  );
  ($($char:expr),* $(,)?) => {
    {
      $crate::representation::seq::Seq::from_iter([$($char),*].map(|c| c as u8).into_iter())
    }
  };
}
