/// Taken with modifications from
/// https://github.com/frangio/concat.rs/blob/d416d3b3c03ba18c0541b7fa63e6e89f3c43e0fe/src/lib.rs
///
/// Credits: Francisco (@frangio)
///
/// Provides the Concat reader adaptor, which wraps around an iterator of readers and exposes its
/// items' contents sequentially. Thus, the contents read from a Concat instance will be the
/// concatenation of the items' contents.
use std::cmp::min;
use std::io::{Read, Result};

pub struct Concat<I>
where
  I: Iterator,
  <I as Iterator>::Item: Read,
{
  iter: I,
  curr: Option<<I as Iterator>::Item>,
  delimiter: Option<Vec<u8>>,
}

impl<I> Concat<I>
where
  I: Iterator,
  <I as Iterator>::Item: Read,
{
  /// Concatenate readers into a single reader
  pub fn from(iter: I) -> Concat<I> {
    Self::with_delimiter(iter, None)
  }

  /// Concatenate readers into a single reader, alternating them with the provided delimiter,
  /// i.e. inserting this sequence of bytes between adjacent readers.
  ///
  /// The idea is that if readers (e.g. files) don't contain trailing newline, it will make line parsers behave
  /// incorrectly on reader boundaries. So one can add newline delimiter between readers, to compensate for that.
  ///
  /// TODO: It might be better to extract the functionality of appending a trailing newline to another iterator class,
  /// and then wrap each reader individually, instead of doing it conditionally in the concatenation iterator like it
  /// is done currently.
  ///
  /// TODO: Ideally, if a reader contains trailing newline(s), we'd like to avoid appending another newline.
  pub fn with_delimiter(mut iter: I, delimiter: Option<Vec<u8>>) -> Concat<I> {
    let curr = iter.next();
    Concat { iter, curr, delimiter }
  }
}

impl<I> Concat<I>
where
  I: Iterator,
  <I as Iterator>::Item: Read,
{
  /// Returns a reference to the item last read, or None if the iterator has been exhausted.
  ///
  /// This is useful for error handling and reporting: if a read operation fails, the reference
  /// returned will point to the item which caused the error.
  #[inline]
  pub const fn current(&self) -> Option<&<I as Iterator>::Item> {
    self.curr.as_ref()
  }
}

impl<I> Read for Concat<I>
where
  I: Iterator,
  <I as Iterator>::Item: Read,
{
  fn read(&mut self, buf: &mut [u8]) -> Result<usize> {
    let n = match self.curr {
      None => 0,
      Some(ref mut r) => r.read(buf)?,
    };

    if n > 0 || buf.is_empty() || self.curr.is_none() {
      Ok(n)
    } else {
      // The current reader reached the end so we have to advance the iterator and try again.
      self.curr = self.iter.next();

      // Before moving to the next reader, insert delimiter, if requested
      let n_bytes_inserted = if let Some(delimiter) = &self.delimiter {
        let n_bytes_inserted = min(delimiter.len(), buf.len());
        buf[..n_bytes_inserted].copy_from_slice(delimiter);
        n_bytes_inserted
      } else {
        0
      };

      let n_bytes_read = self.read(&mut buf[n_bytes_inserted..])?;
      Ok(n_bytes_read + n_bytes_inserted)
    }
  }
}
