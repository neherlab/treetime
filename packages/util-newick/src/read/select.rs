use crate::dialect::NewickDialect;
use crate::grammar::{Rule, parse, rule_name, start_rule};
use crate::read::builder::build_graph;
use crate::read::context::MapContext;
use crate::read::error::{DialectAttempt, Location, NewickError, NewickErrorKind, TextIndex};
use crate::read::options::{NewickReadOptions, NewickTree, ReadMode};
use pest::error::{Error, ErrorVariant, InputLocation};

pub(crate) fn read_tree_text(
  text: &str,
  base: Location,
  options: &NewickReadOptions,
) -> Result<NewickTree, NewickError> {
  if options.dialects.is_empty() {
    return Err(NewickError::new(
      NewickErrorKind::Options,
      base,
      "No Newick dialect is selected",
    ));
  }
  let index = TextIndex::new(text, base);
  let modes: &[ReadMode] = match options.mode {
    ReadMode::Strict => &[ReadMode::Strict],
    ReadMode::Tolerant => &[ReadMode::Strict, ReadMode::Tolerant],
  };
  let mut attempts: Vec<DialectAttempt> = Vec::new();
  for &mode in modes {
    for &dialect in &options.dialects {
      match read_in_dialect(text, &index, dialect, mode, options) {
        Ok(tree) => return Ok(tree),
        Err(error) => match attempts.iter_mut().find(|attempt| attempt.dialect == dialect) {
          Some(attempt) => {
            attempt.mode = mode;
            attempt.error = error;
          },
          None => attempts.push(DialectAttempt { dialect, mode, error }),
        },
      }
    }
  }
  match <[DialectAttempt; 1]>::try_from(attempts) {
    Ok([attempt]) => Err(attempt.error),
    Err(attempts) => Err(NewickError::new(
      NewickErrorKind::NoDialect(attempts),
      base,
      "No selected Newick dialect reads the tree:",
    )),
  }
}

pub(crate) fn syntax_error(error: &Error<Rule>, index: &TextIndex<'_>) -> NewickError {
  let offset = match error.location {
    InputLocation::Pos(offset) | InputLocation::Span((offset, _)) => offset,
  };
  let message = match &error.variant {
    ErrorVariant::ParsingError { positives, negatives } => {
      let expected = rule_names(positives);
      let unexpected = rule_names(negatives);
      match (expected.is_empty(), unexpected.is_empty()) {
        (false, true) => format!("expected {expected}"),
        (true, false) => format!("unexpected {unexpected}"),
        (false, false) => format!("expected {expected}, unexpected {unexpected}"),
        (true, true) => "unexpected input".to_owned(),
      }
    },
    ErrorVariant::CustomError { message } => message.clone(),
  };
  NewickError::new(NewickErrorKind::Syntax, index.locate(offset), message)
}

fn read_in_dialect(
  text: &str,
  index: &TextIndex<'_>,
  dialect: NewickDialect,
  mode: ReadMode,
  options: &NewickReadOptions,
) -> Result<NewickTree, NewickError> {
  let tree = parse(start_rule(dialect, mode), text).map_err(|error| syntax_error(&error, index))?;
  let mut context = MapContext::new(options, mode, index);
  let tokens = tree.flat_map(|start| start.into_inner());
  let graph = build_graph(tokens, &mut context, index.locate(0))?;
  Ok(NewickTree {
    graph,
    dialect,
    warnings: context.warnings,
  })
}

fn rule_names(rules: &[Rule]) -> String {
  let mut names: Vec<&str> = Vec::new();
  for name in rules.iter().map(|&rule| rule_name(rule)) {
    if !names.contains(&name) {
      names.push(name);
    }
  }
  names.join(", ")
}
