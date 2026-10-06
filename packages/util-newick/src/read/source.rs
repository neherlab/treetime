use crate::read::error::{Location, NewickError, NewickErrorKind};
use std::io::{self, Read};
use std::str;

const FIRST_READ_SIZE: usize = 64 * 1024;

pub(crate) struct TextSource<R> {
  reader: R,
  buffer: Vec<u8>,
  start: usize,
  location: Location,
  read_size: usize,
  eof: bool,
}

impl<R: Read> TextSource<R> {
  pub(crate) fn new(reader: R) -> Self {
    Self {
      reader,
      buffer: Vec::new(),
      start: 0,
      location: Location::START,
      read_size: FIRST_READ_SIZE,
      eof: false,
    }
  }

  pub(crate) fn location(&self) -> Location {
    self.location
  }

  pub(crate) fn window(&self) -> Window<'_> {
    let available = self.buffer.get(self.start..).unwrap_or_default();
    match str::from_utf8(available) {
      Ok(text) => Window {
        text,
        at_end: self.eof,
        is_invalid: false,
      },
      Err(error) => {
        let valid = available.get(..error.valid_up_to()).unwrap_or_default();
        let text = str::from_utf8(valid).unwrap_or_default();
        let waits_for_more = error.error_len().is_none() && !self.eof;
        Window {
          text,
          at_end: false,
          is_invalid: !waits_for_more,
        }
      },
    }
  }

  pub(crate) fn advance(&mut self, consumed_len: usize, location: Location) {
    self.start += consumed_len;
    self.location = location;
  }

  pub(crate) fn fill(&mut self) -> Result<(), NewickError> {
    let window = self.window();
    if window.is_invalid {
      let location = self.location.advanced_by(window.text);
      return Err(NewickError::new(
        NewickErrorKind::InvalidUtf8,
        location,
        "The input is not valid UTF-8",
      ));
    }
    self.read_more().map_err(|error| {
      NewickError::new(
        NewickErrorKind::Io,
        self.location,
        format!("When reading the input: {error}"),
      )
    })
  }

  fn read_more(&mut self) -> io::Result<()> {
    if self.start > 0 && self.start * 2 >= self.buffer.len() {
      self.buffer.drain(..self.start);
      self.start = 0;
    }
    let limit = u64::try_from(self.read_size).unwrap_or(u64::MAX);
    let read = self.reader.by_ref().take(limit).read_to_end(&mut self.buffer)?;
    self.eof = read == 0;
    self.read_size = self.read_size.saturating_mul(2);
    Ok(())
  }
}

pub(crate) struct Window<'b> {
  pub(crate) text: &'b str,
  pub(crate) at_end: bool,
  pub(crate) is_invalid: bool,
}
