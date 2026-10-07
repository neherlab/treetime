use crate::dialect::NewickAnnotations;
use crate::grammar::{Rule, matches, parse};
use crate::read::error::{Location, NewickError, NewickErrorKind};
use crate::read::options::{NewickReadOptions, NewickTree};
use crate::read::source::{Scan, TextSource, scan_to_semicolon};
use crate::read::tree_text::read_tree_text;
use std::io::Read;
use std::str;

pub fn newick_from_str(input: &str, options: &NewickReadOptions) -> Result<NewickTree, NewickError> {
  if matches(Rule::trivia_only, input) {
    return Err(no_tree());
  }
  read_tree_text(input, Location::START, options).map_err(|error| {
    let mut trees = newick_trees(input.as_bytes(), options.clone());
    match (trees.next(), trees.next()) {
      (Some(Err(scan_error)), _) if scan_error.kind == NewickErrorKind::Incomplete => scan_error,
      (Some(_), Some(_)) => NewickError::new(
        NewickErrorKind::MultipleTrees,
        trees.last_tree_start,
        "The input contains more than one tree; read it tree by tree with newick_trees()",
      ),
      _ => error,
    }
  })
}

pub fn newick_from_reader(mut reader: impl Read, options: &NewickReadOptions) -> Result<NewickTree, NewickError> {
  let mut bytes = Vec::new();
  reader.read_to_end(&mut bytes).map_err(|error| {
    NewickError::new(
      NewickErrorKind::Io,
      Location::START,
      format!("When reading the input: {error}"),
    )
  })?;
  match str::from_utf8(&bytes) {
    Ok(input) => newick_from_str(input, options),
    Err(error) => {
      let valid = bytes.get(..error.valid_up_to()).unwrap_or_default();
      let location = Location::START.advanced_by(str::from_utf8(valid).unwrap_or_default());
      Err(NewickError::new(
        NewickErrorKind::InvalidUtf8,
        location,
        "The input is not valid UTF-8",
      ))
    },
  }
}

pub fn newick_trees<R: Read>(reader: R, options: NewickReadOptions) -> NewickTrees<R> {
  NewickTrees {
    source: TextSource::new(reader),
    options,
    scanned: 0,
    last_tree_start: Location::START,
    done: false,
  }
}

pub struct NewickTrees<R> {
  source: TextSource<R>,
  options: NewickReadOptions,
  scanned: usize,
  last_tree_start: Location,
  done: bool,
}

impl<R: Read> NewickTrees<R> {
  fn step(&mut self) -> Option<Result<NewickTree, NewickError>> {
    loop {
      let window = self.source.window();
      let location = self.source.location();
      self.last_tree_start = location;
      match scan_to_semicolon(
        window.text,
        window.at_end,
        self.options.mode,
        &mut self.scanned,
        |text, at_end| scan_prefix(text, at_end, self.options.dialect.annotations),
        |text| matches(Rule::trivia_only, text),
      ) {
        Scan::Found(slice) => {
          let result = read_tree_text(slice, location, &self.options);
          let (consumed_len, next) = (slice.len(), location.advanced_by(slice));
          self.source.advance(consumed_len, next);
          return Some(result);
        },
        Scan::Unterminated(rest) => {
          let at = location.advanced_by(rest.trim_end());
          return Some(Err(NewickError::new(
            NewickErrorKind::Incomplete,
            at,
            "The tree does not end with ';'",
          )));
        },
        Scan::Trivia => return None,
        Scan::NeedMore => {},
      }
      if let Err(error) = self.source.fill() {
        return Some(Err(error));
      }
    }
  }
}

impl<R: Read> Iterator for NewickTrees<R> {
  type Item = Result<NewickTree, NewickError>;

  fn next(&mut self) -> Option<Self::Item> {
    if self.done {
      return None;
    }
    let item = self.step();
    if !matches!(item, Some(Ok(_))) {
      self.done = true;
    }
    item
  }
}

fn scan_prefix(text: &str, at_end: bool, annotations: NewickAnnotations) -> usize {
  let rule = match (annotations.reserves_annotations(), at_end) {
    (true, false) => Rule::tree_scan,
    (true, true) => Rule::tree_scan_final,
    (false, false) => Rule::tree_scan_plain,
    (false, true) => Rule::tree_scan_plain_final,
  };
  parse(rule, text)
    .ok()
    .and_then(|mut pairs| pairs.next())
    .map_or(0, |prefix| prefix.as_str().len())
}

fn no_tree() -> NewickError {
  NewickError::new(NewickErrorKind::Syntax, Location::START, "The input contains no tree")
}
