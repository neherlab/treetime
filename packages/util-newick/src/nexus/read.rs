use crate::dialect::NewickDialect;
use crate::grammar;
use crate::grammar::comment_rule;
use crate::model::comment::NewickComment;
use crate::model::graph::NewickGraph;
use crate::nexus::grammar::{Rule, matches, parse, rule_name};
use crate::nexus::types::{NexusCommand, NexusFile, NexusTree};
use crate::read::comments::{read_comment, unquote};
use crate::read::context::MapContext;
use crate::read::error::{Location, NewickError, NewickErrorKind, NewickWarning, TextIndex};
use crate::read::options::{NewickReadOptions, ReadMode};
use crate::read::select::{read_tree_text, syntax_error};
use crate::read::source::{Scan, TextSource, scan_to_semicolon};
use pest::error::{Error, ErrorVariant, InputLocation};
use pest::iterators::Pair;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Read;

pub fn is_nexus(input: &str) -> bool {
  matches(Rule::nexus_header, input)
}

pub fn nexus_from_str(input: &str, options: &NewickReadOptions) -> Result<NexusFile, NewickError> {
  nexus_from_reader(input.as_bytes(), options)
}

pub fn nexus_from_reader(reader: impl Read, options: &NewickReadOptions) -> Result<NexusFile, NewickError> {
  let mut reading = nexus_trees(reader, options.clone());
  let trees = reading.by_ref().collect::<Result<Vec<_>, _>>()?;
  Ok(NexusFile {
    trees,
    skipped: reading.skipped,
    warnings: reading.warnings,
  })
}

pub fn nexus_trees<R: Read>(reader: R, options: NewickReadOptions) -> NexusTrees<R> {
  NexusTrees {
    source: TextSource::new(reader),
    options,
    state: BlockState::default(),
    scanned: 0,
    has_header: false,
    done: false,
    skipped: Vec::new(),
    warnings: Vec::new(),
  }
}

pub struct NexusTrees<R> {
  source: TextSource<R>,
  options: NewickReadOptions,
  state: BlockState,
  scanned: usize,
  has_header: bool,
  done: bool,
  skipped: Vec<NexusCommand>,
  warnings: Vec<NewickWarning>,
}

impl<R: Read> NexusTrees<R> {
  pub fn skipped(&self) -> &[NexusCommand] {
    &self.skipped
  }

  pub fn warnings(&self) -> &[NewickWarning] {
    &self.warnings
  }

  fn step(&mut self) -> Option<Result<NexusTree, NewickError>> {
    loop {
      let window = self.source.window();
      let location = self.source.location();
      if !self.has_header {
        match read_header(window.text, window.at_end, location) {
          Header::Found(len) => {
            let next = location.advanced_by(window.text.get(..len).unwrap_or_default());
            self.source.advance(len, next);
            self.has_header = true;
            continue;
          },
          Header::Missing(error) => return Some(Err(error)),
          Header::NeedMore => {},
        }
      } else {
        match scan_to_semicolon(
          window.text,
          window.at_end,
          self.options.mode,
          &mut self.scanned,
          scan_prefix,
          |text| matches(Rule::trivia_only, text),
        ) {
          Scan::Found(slice) => {
            let (len, next) = (slice.len(), location.advanced_by(slice));
            let mut reader = CommandReader {
              options: &self.options,
              state: &mut self.state,
              skipped: &mut self.skipped,
              warnings: &mut self.warnings,
            };
            let outcome = reader.read(slice, location);
            self.source.advance(len, next);
            match outcome {
              Ok(Some(tree)) => return Some(Ok(tree)),
              Ok(None) => continue,
              Err(error) => return Some(Err(error)),
            }
          },
          Scan::Unterminated(rest) => {
            let at = location.advanced_by(rest.trim_end());
            return Some(Err(NewickError::new(
              NewickErrorKind::Incomplete,
              at,
              "The command does not end with ';'",
            )));
          },
          Scan::Trivia => return self.finish(location.advanced_by(window.text)).err().map(Err),
          Scan::NeedMore => {},
        }
      }
      if let Err(error) = self.source.fill() {
        return Some(Err(error));
      }
    }
  }

  fn finish(&mut self, location: Location) -> Result<(), NewickError> {
    let Some(block) = &self.state.block else {
      return Ok(());
    };
    let message = format!("The {} block never ends", block.name);
    tolerate(self.options.mode, &mut self.warnings, location, message)
  }
}

impl<R: Read> Iterator for NexusTrees<R> {
  type Item = Result<NexusTree, NewickError>;

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

#[derive(Default)]
struct BlockState {
  block: Option<Block>,
  translate: BTreeMap<String, String>,
  rooted: Option<bool>,
  taxa: Option<Taxa>,
  ntax: Option<usize>,
}

struct Block {
  name: String,
  kind: BlockKind,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BlockKind {
  Trees,
  Taxa,
  Other,
}

struct Taxa {
  labels: Vec<String>,
  known: BTreeSet<String>,
}

enum Header {
  Found(usize),
  Missing(NewickError),
  NeedMore,
}

fn read_header(text: &str, at_end: bool, location: Location) -> Header {
  if let Ok(mut pairs) = parse(Rule::nexus_header, text)
    && let Some(header) = pairs.next()
  {
    return Header::Found(header.as_span().end());
  }
  if !at_end && text.trim_start().len() < "#NEXUS ".len() {
    return Header::NeedMore;
  }
  Header::Missing(NewickError::new(
    NewickErrorKind::Nexus,
    location,
    "The input does not start with a '#NEXUS' header",
  ))
}

fn scan_prefix(text: &str, at_end: bool) -> usize {
  let rule = if at_end {
    Rule::command_scan_final
  } else {
    Rule::command_scan
  };
  parse(rule, text)
    .ok()
    .and_then(|mut pairs| pairs.next())
    .map_or(0, |prefix| prefix.as_str().len())
}

struct CommandReader<'r> {
  options: &'r NewickReadOptions,
  state: &'r mut BlockState,
  skipped: &'r mut Vec<NexusCommand>,
  warnings: &'r mut Vec<NewickWarning>,
}

impl CommandReader<'_> {
  fn read(&mut self, slice: &str, location: Location) -> Result<Option<NexusTree>, NewickError> {
    let index = TextIndex::new(slice, location);
    let rule = match self.options.mode {
      ReadMode::Strict => Rule::command_strict,
      ReadMode::Tolerant => Rule::command_tolerant,
    };
    let pairs = parse(rule, slice).map_err(|error| nexus_syntax_error(&error, &index))?;
    let Some(command) = pairs
      .flat_map(|pair| pair.into_inner())
      .find(|pair| !matches!(pair.as_rule(), Rule::scan_comment | Rule::terminator | Rule::EOI))
    else {
      return Ok(None);
    };
    let at = index.locate(command.as_span().start());
    let kind = self.state.block.as_ref().map(|block| block.kind);
    match (command.as_rule(), kind) {
      (Rule::begin_command, _) => self.begin(&command, at),
      (Rule::end_command, _) => self.end(at),
      (Rule::tree_strict | Rule::tree_tolerant, Some(BlockKind::Trees)) => {
        let newick_start = command
          .clone()
          .into_inner()
          .find(|part| part.as_rule() == Rule::newick_text)
          .map_or(slice.len(), |part| part.as_span().start());
        let newick = slice.get(newick_start..).unwrap_or_default();
        return self.read_tree(command, newick, &index, newick_start, at).map(Some);
      },
      (Rule::translate_strict | Rule::translate_tolerant, Some(BlockKind::Trees)) => self.translate(command, at),
      (Rule::properties_command, Some(BlockKind::Trees)) => {
        self.properties(command);
        Ok(())
      },
      (Rule::taxlabels_command, Some(BlockKind::Taxa)) => self.taxlabels(command, at),
      (Rule::dimensions_command, Some(BlockKind::Taxa)) => {
        self.dimensions(command);
        Ok(())
      },
      (_, kind) => {
        if kind.is_none() {
          tolerate(
            self.options.mode,
            self.warnings,
            at,
            "A command appears outside of a block",
          )?;
        }
        self.skip(&command, at);
        Ok(())
      },
    }?;
    Ok(None)
  }

  fn begin(&mut self, command: &Pair<'_, Rule>, at: Location) -> Result<(), NewickError> {
    if let Some(block) = &self.state.block {
      let message = format!("A block begins before the {} block ends", block.name);
      tolerate(self.options.mode, self.warnings, at, message)?;
    }
    let name = command
      .clone()
      .into_inner()
      .find(|part| part.as_rule() == Rule::word)
      .map(|word| self.word(&word))
      .unwrap_or_default();
    let kind = match name.to_ascii_lowercase().as_str() {
      "trees" => BlockKind::Trees,
      "taxa" => BlockKind::Taxa,
      _ => BlockKind::Other,
    };
    self.state.block = Some(Block { name, kind });
    self.state.translate.clear();
    self.state.rooted = None;
    Ok(())
  }

  fn end(&mut self, at: Location) -> Result<(), NewickError> {
    if self.state.block.take().is_none() {
      tolerate(self.options.mode, self.warnings, at, "'End' appears outside of a block")?;
    }
    Ok(())
  }

  fn skip(&mut self, command: &Pair<'_, Rule>, at: Location) {
    let keyword = command
      .clone()
      .into_inner()
      .next()
      .map_or_else(|| command.as_str().to_owned(), |first| first.as_str().to_owned());
    self.skipped.push(NexusCommand {
      block: self.state.block.as_ref().map(|block| block.name.clone()),
      command: keyword,
      line: at.line,
    });
  }

  fn word(&self, word: &Pair<'_, Rule>) -> String {
    let inner = word.clone().into_inner().next();
    let token = inner.as_ref().unwrap_or(word);
    match token.as_rule() {
      Rule::quoted_label => unquote(token.as_str(), '\''),
      _ if self.options.underscores_as_spaces => token.as_str().replace('_', " "),
      _ => token.as_str().to_owned(),
    }
  }

  fn translate(&mut self, command: Pair<'_, Rule>, at: Location) -> Result<(), NewickError> {
    for pair in command
      .into_inner()
      .filter(|part| part.as_rule() != Rule::translate_keyword)
    {
      let words: Vec<String> = pair
        .into_inner()
        .filter(|part| !matches!(part.as_rule(), Rule::scan_comment | Rule::separator))
        .map(|word| self.word(&word))
        .collect();
      let [key, value] = words.as_slice() else {
        continue;
      };
      if self.state.translate.insert(key.clone(), value.clone()).is_some() {
        let message = format!("The 'Translate' command lists the key {key:?} more than once");
        tolerate(self.options.mode, self.warnings, at, message)?;
      }
    }
    Ok(())
  }

  fn properties(&mut self, command: Pair<'_, Rule>) {
    for property in command.into_inner().filter(|part| part.as_rule() == Rule::property) {
      let words: Vec<String> = property
        .into_inner()
        .filter(|part| part.as_rule() == Rule::word)
        .map(|word| self.word(&word))
        .collect();
      if let [key, value] = words.as_slice()
        && key.eq_ignore_ascii_case("rooted")
      {
        self.state.rooted = Some(matches!(value.to_ascii_lowercase().as_str(), "yes" | "true"));
      }
    }
  }

  fn taxlabels(&mut self, command: Pair<'_, Rule>, at: Location) -> Result<(), NewickError> {
    let labels: Vec<String> = command
      .into_inner()
      .filter(|part| part.as_rule() == Rule::word)
      .map(|word| self.word(&word))
      .collect();
    if let Some(ntax) = self.state.ntax
      && ntax != labels.len()
    {
      let message = format!(
        "The 'Dimensions' command declares {ntax} taxa, but 'TaxLabels' lists {}",
        labels.len()
      );
      tolerate(self.options.mode, self.warnings, at, message)?;
    }
    let known = labels.iter().cloned().collect();
    self.state.taxa = Some(Taxa { labels, known });
    Ok(())
  }

  fn dimensions(&mut self, command: Pair<'_, Rule>) {
    let words: Vec<String> = command
      .into_inner()
      .filter(|part| part.as_rule() == Rule::word)
      .map(|word| self.word(&word))
      .collect();
    for pair in words.chunks(2) {
      if let [key, value] = pair
        && key.eq_ignore_ascii_case("ntax")
        && let Ok(ntax) = value.parse::<usize>()
      {
        self.state.ntax = Some(ntax);
      }
    }
  }

  fn read_tree(
    &self,
    command: Pair<'_, Rule>,
    newick: &str,
    index: &TextIndex<'_>,
    newick_start: usize,
    at: Location,
  ) -> Result<NexusTree, NewickError> {
    let mut name = String::new();
    let mut comment_texts = Vec::new();
    for part in command.into_inner() {
      match part.as_rule() {
        Rule::word => name = self.word(&part),
        Rule::quoted_label => name = unquote(part.as_str(), '\''),
        Rule::loose_text => name = part.as_str().to_owned(),
        Rule::scan_comment if !name.is_empty() => {
          comment_texts.push((part.as_str(), index.locate(part.as_span().start())));
        },
        _ => {},
      }
    }
    let mut tree = read_tree_text(newick, index.locate(newick_start), self.options)?;
    let mut comments = Vec::new();
    for (text, location) in comment_texts {
      let comment = tree_comment(text, location, &mut tree.warnings, self.options, tree.dialect)?;
      comments.push(comment);
    }
    self.resolve_leaf_names(&mut tree.graph, &mut tree.warnings, at)?;
    if tree.graph.rooted().is_none() {
      tree.graph.set_rooted(self.state.rooted);
    }
    Ok(NexusTree { name, tree, comments })
  }

  fn resolve_leaf_names(
    &self,
    graph: &mut NewickGraph,
    warnings: &mut Vec<NewickWarning>,
    at: Location,
  ) -> Result<(), NewickError> {
    for node in 0..graph.node_count() {
      if !graph.is_leaf(node) {
        continue;
      }
      let Some(name) = graph.node(node).name() else {
        continue;
      };
      let resolved = match (self.state.translate.get(name), &self.state.taxa) {
        (Some(translated), _) => Some(translated.clone()),
        (None, Some(taxa)) if taxa.known.contains(name) => None,
        (None, Some(taxa)) => {
          if let Some(label) = name
            .parse::<usize>()
            .ok()
            .and_then(|number| taxa.labels.get(number.checked_sub(1)?))
          {
            Some(label.clone())
          } else {
            let message = format!("The tree label {name:?} is neither a taxon label nor a 'Translate' key");
            tolerate(self.options.mode, warnings, at, message)?;
            None
          }
        },
        (None, None) => None,
      };
      if let Some(resolved) = resolved {
        graph.node_mut(node).set_name(Some(resolved));
      }
    }
    Ok(())
  }
}

fn tree_comment(
  text: &str,
  at: Location,
  warnings: &mut Vec<NewickWarning>,
  options: &NewickReadOptions,
  dialect: NewickDialect,
) -> Result<NewickComment, NewickError> {
  let index = TextIndex::new(text, at);
  let mut pairs = grammar::parse(comment_rule(dialect), text).map_err(|error| syntax_error(&error, &index))?;
  let mut context = MapContext::new(options, options.mode, &index);
  let comment = pairs
    .next()
    .and_then(|only| only.into_inner().next())
    .map(|comment| read_comment(comment, &mut context))
    .transpose()?
    .unwrap_or_else(|| NewickComment::Plain(String::new()));
  warnings.append(&mut context.warnings);
  Ok(comment)
}

fn tolerate(
  mode: ReadMode,
  warnings: &mut Vec<NewickWarning>,
  location: Location,
  message: impl Into<String>,
) -> Result<(), NewickError> {
  match mode {
    ReadMode::Strict => Err(NewickError::new(NewickErrorKind::Nexus, location, message)),
    ReadMode::Tolerant => {
      warnings.push(NewickWarning::new(location, message));
      Ok(())
    },
  }
}

fn nexus_syntax_error(error: &Error<Rule>, index: &TextIndex<'_>) -> NewickError {
  let offset = match error.location {
    InputLocation::Pos(offset) | InputLocation::Span((offset, _)) => offset,
  };
  let message = match &error.variant {
    ErrorVariant::ParsingError { positives, .. } => {
      let mut names: Vec<&str> = Vec::new();
      for name in positives.iter().map(|&rule| rule_name(rule)) {
        if !names.contains(&name) {
          names.push(name);
        }
      }
      format!("expected {}", names.join(", "))
    },
    ErrorVariant::CustomError { message } => message.clone(),
  };
  NewickError::new(NewickErrorKind::Syntax, index.locate(offset), message)
}
