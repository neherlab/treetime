use crate::grammar::{Rule, parse, rule_name, start_rule};
use crate::read::builder::build_graph;
use crate::read::context::MapContext;
use crate::read::error::{Location, NewickError, NewickErrorKind, TextIndex};
use crate::read::options::{NewickReadOptions, NewickTree};
use pest::error::{Error, ErrorVariant, InputLocation};

pub(crate) fn read_tree_text(
  text: &str,
  base: Location,
  options: &NewickReadOptions,
) -> Result<NewickTree, NewickError> {
  let index = TextIndex::new(text, base);
  let tree = parse(start_rule(options.dialect, options.mode), text).map_err(|error| syntax_error(&error, &index))?;
  let mut context = MapContext::new(options, &index);
  let tokens = tree.flat_map(|start| start.into_inner());
  let graph = build_graph(tokens, &mut context, index.locate(0))?;
  Ok(NewickTree {
    graph,
    warnings: context.warnings,
  })
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

fn rule_names(rules: &[Rule]) -> String {
  let mut names: Vec<&str> = Vec::new();
  for name in rules.iter().map(|&rule| rule_name(rule)) {
    if !names.contains(&name) {
      names.push(name);
    }
  }
  names.join(", ")
}
