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
  pub(crate) fn with_delimiter(mut iter: I, delimiter: Option<Vec<u8>>) -> Concat<I> {
    let curr = iter.next();
    Concat { iter, curr, delimiter }
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
      self.curr = self.iter.next();

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
