use std::error::Error;
use std::fmt;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NewickError {
  pub kind: NewickErrorKind,
  pub offset: usize,
  pub line: usize,
  pub column: usize,
  pub message: String,
}

impl NewickError {
  pub(crate) fn new(kind: NewickErrorKind, location: Location, message: impl Into<String>) -> Self {
    Self {
      kind,
      offset: location.offset,
      line: location.line,
      column: location.column,
      message: message.into(),
    }
  }

  pub fn location(&self) -> Location {
    Location {
      offset: self.offset,
      line: self.line,
      column: self.column,
    }
  }
}

impl fmt::Display for NewickError {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "line {}, column {}: {}", self.line, self.column, self.message)
  }
}

impl Error for NewickError {}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NewickErrorKind {
  Syntax,
  Structure,
  Annotation,
  Nexus,
  InvalidUtf8,
  MultipleTrees,
  Incomplete,
  Io,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NewickWarning {
  pub offset: usize,
  pub line: usize,
  pub column: usize,
  pub message: String,
}

impl NewickWarning {
  pub(crate) fn new(location: Location, message: impl Into<String>) -> Self {
    Self {
      offset: location.offset,
      line: location.line,
      column: location.column,
      message: message.into(),
    }
  }
}

impl fmt::Display for NewickWarning {
  fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    write!(f, "line {}, column {}: {}", self.line, self.column, self.message)
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Location {
  pub offset: usize,
  pub line: usize,
  pub column: usize,
}

impl Location {
  pub const START: Self = Self {
    offset: 0,
    line: 1,
    column: 1,
  };

  pub(crate) fn advanced_by(self, text: &str) -> Self {
    let offset = self.offset + text.len();
    match text.rfind('\n') {
      Some(last_newline) => Self {
        offset,
        line: self.line + text.matches('\n').count(),
        column: text.get(last_newline + 1..).map_or(0, |rest| rest.chars().count()) + 1,
      },
      None => Self {
        offset,
        line: self.line,
        column: self.column + text.chars().count(),
      },
    }
  }
}

pub(crate) struct TextIndex<'i> {
  text: &'i str,
  base: Location,
  line_starts: Vec<usize>,
}

impl<'i> TextIndex<'i> {
  pub(crate) fn new(text: &'i str, base: Location) -> Self {
    let line_starts = std::iter::once(0)
      .chain(text.match_indices('\n').map(|(idx, _)| idx + 1))
      .collect();
    Self {
      text,
      base,
      line_starts,
    }
  }

  pub(crate) fn locate(&self, offset: usize) -> Location {
    let line = self
      .line_starts
      .partition_point(|&start| start <= offset)
      .saturating_sub(1);
    let line_start = self.line_starts.get(line).copied().unwrap_or(0);
    let column = self
      .text
      .get(line_start..offset)
      .map_or(0, |prefix| prefix.chars().count());
    if line == 0 {
      Location {
        offset: self.base.offset + offset,
        line: self.base.line,
        column: self.base.column + column,
      }
    } else {
      Location {
        offset: self.base.offset + offset,
        line: self.base.line + line,
        column: column + 1,
      }
    }
  }
}
