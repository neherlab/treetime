use crate::grammar::Rule;
use crate::read::error::{Location, NewickError, NewickErrorKind, NewickWarning, TextIndex};
use crate::read::options::{NewickReadOptions, ReadMode};
use pest::iterators::Pair;

pub(crate) struct MapContext<'o, 'i> {
  pub(crate) options: &'o NewickReadOptions,
  pub(crate) index: &'o TextIndex<'i>,
  pub(crate) warnings: Vec<NewickWarning>,
}

impl<'o, 'i> MapContext<'o, 'i> {
  pub(crate) fn new(options: &'o NewickReadOptions, index: &'o TextIndex<'i>) -> Self {
    Self {
      options,
      index,
      warnings: Vec::new(),
    }
  }

  pub(crate) fn locate(&self, pair: &Pair<'_, Rule>) -> Location {
    self.index.locate(pair.as_span().start())
  }

  pub(crate) fn error(&self, kind: NewickErrorKind, pair: &Pair<'_, Rule>, message: impl Into<String>) -> NewickError {
    NewickError::new(kind, self.locate(pair), message)
  }

  pub(crate) fn tolerate(
    &mut self,
    kind: NewickErrorKind,
    pair: &Pair<'_, Rule>,
    message: impl Into<String>,
  ) -> Result<(), NewickError> {
    let location = self.locate(pair);
    self.tolerate_at(kind, location, message)
  }

  pub(crate) fn tolerate_at(
    &mut self,
    kind: NewickErrorKind,
    location: Location,
    message: impl Into<String>,
  ) -> Result<(), NewickError> {
    match self.options.mode {
      ReadMode::Strict => Err(NewickError::new(kind, location, message)),
      ReadMode::Tolerant => {
        self.warnings.push(NewickWarning::new(location, message));
        Ok(())
      },
    }
  }
}
