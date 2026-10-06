use crate::grammar::{Rule, matches, parse};
use crate::read::error::{Location, NewickError, NewickErrorKind};
use crate::read::options::{NewickReadOptions, NewickTree, ReadMode};
use crate::read::select::read_tree_text;
use crate::read::source::TextSource;
use std::io::Read;

pub fn newick_from_str(input: &str, options: &NewickReadOptions) -> Result<NewickTree, NewickError> {
  newick_from_reader(input.as_bytes(), options)
}

pub fn newick_from_reader(reader: impl Read, options: &NewickReadOptions) -> Result<NewickTree, NewickError> {
  let mut trees = newick_trees(reader, options.clone());
  let tree = trees.next().ok_or_else(no_tree)??;
  let location = trees.source.location();
  match trees.next() {
    None => Ok(tree),
    Some(_) => Err(NewickError::new(
      NewickErrorKind::MultipleTrees,
      location,
      "The input contains more than one tree; read it tree by tree with newick_trees()",
    )),
  }
}

pub fn newick_trees<R: Read>(reader: R, options: NewickReadOptions) -> NewickTrees<R> {
  NewickTrees {
    source: TextSource::new(reader),
    options,
    done: false,
  }
}

pub struct NewickTrees<R> {
  source: TextSource<R>,
  options: NewickReadOptions,
  done: bool,
}

impl<R: Read> NewickTrees<R> {
  fn step(&mut self) -> Option<Result<NewickTree, NewickError>> {
    loop {
      let window = self.source.window();
      let location = self.source.location();
      match scan_tree(window.text, window.at_end, self.options.mode) {
        TreeScan::Tree(slice) => {
          let result = read_tree_text(slice, location, &self.options);
          let (consumed_len, next) = (slice.len(), location.advanced_by(slice));
          self.source.advance(consumed_len, next);
          return Some(result);
        },
        TreeScan::Unterminated(rest) => {
          let at = location.advanced_by(rest.trim_end());
          return Some(Err(NewickError::new(
            NewickErrorKind::Syntax,
            at,
            "The tree does not end with ';'",
          )));
        },
        TreeScan::Done => return None,
        TreeScan::NeedMore => {},
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

enum TreeScan<'t> {
  Tree(&'t str),
  Unterminated(&'t str),
  Done,
  NeedMore,
}

fn scan_tree(text: &str, at_end: bool, mode: ReadMode) -> TreeScan<'_> {
  let rule = if at_end {
    Rule::tree_extent_final
  } else {
    Rule::tree_extent
  };
  if let Ok(mut pairs) = parse(rule, text)
    && let Some(extent) = pairs.next()
  {
    return TreeScan::Tree(extent.as_str());
  }
  match (at_end, mode) {
    (false, _) => TreeScan::NeedMore,
    (true, _) if matches(Rule::trivia_only, text) => TreeScan::Done,
    (true, ReadMode::Strict) => TreeScan::Unterminated(text),
    (true, ReadMode::Tolerant) => TreeScan::Tree(text),
  }
}

fn no_tree() -> NewickError {
  NewickError::new(NewickErrorKind::Syntax, Location::START, "The input contains no tree")
}
