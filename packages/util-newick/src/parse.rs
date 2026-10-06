use crate::annotation::classify_comment;
use crate::types::{NewickEdgeData, NewickGraph, NewickHybrid, NewickLabel, NewickNodeData, NewickReadOptions};
use crate::validate::describe_node;
use eyre::{Report, WrapErr, eyre};
use pest::Parser;
use pest::error::Error as PestError;
use pest::iterators::Pair;
use pest_derive::Parser;
use regex::regex;
use std::collections::BTreeMap;
use std::io::Read;

pub fn newick_from_reader(mut reader: impl Read, options: &NewickReadOptions) -> Result<NewickGraph, Report> {
  let mut input = String::new();
  reader
    .read_to_string(&mut input)
    .wrap_err("When reading Newick input")?;
  newick_from_string(&input, options)
}

pub fn newick_from_string(input: &str, options: &NewickReadOptions) -> Result<NewickGraph, Report> {
  build_graph(input, options).wrap_err("Failed to parse Newick string")
}

pub(crate) fn is_comment_token(text: &str) -> bool {
  NewickParser::parse(Rule::comment_token, text).is_ok()
}

#[derive(Parser)]
#[grammar_inline = r##"
tree           = _{ SOI ~ "\u{FEFF}"? ~ rooting? ~ item* ~ (end ~ comment*)? ~ EOI }
item           = _{ open | close | comma | length | label | comment }
open           =  { "(" }
close          =  { ")" }
comma          =  { "," }
end            =  { ";" }
length         =  { ":" ~ comment* ~ number? }
rooting        =  { "[&" ~ rooting_value ~ "]" }
rooting_value  = @{ ^"r" | ^"u" }
label          = ${ quoted_label ~ hybrid_tag? | unquoted_label }
quoted_label   = @{ "'" ~ ("''" | !"'" ~ ANY)* ~ "'" }
hybrid_tag     = @{ "#" ~ "#"? ~ ASCII_ALPHA* ~ ASCII_DIGIT+ }
unquoted_label = @{ !"'" ~ (!("(" | ")" | "[" | "]" | "," | ";" | ":" | WHITESPACE) ~ ANY)+ }
number         = @{ ("+" | "-")? ~ (ASCII_DIGIT+ ~ ("." ~ ASCII_DIGIT*)? | "." ~ ASCII_DIGIT+) ~ (^"e" ~ ("+" | "-")? ~ ASCII_DIGIT+)? }
comment        = @{ "[" ~ comment_char* ~ "]" }
comment_char   = _{ "\"" ~ (!"\"" ~ ANY)* ~ "\"" | PUSH("[") | "]" ~ DROP | !"]" ~ ANY }
comment_token  = _{ SOI ~ strict_comment ~ EOI }
strict_comment = @{ "[" ~ ("\"" ~ (!"\"" ~ ANY)* ~ "\"" | PUSH("[") | "]" ~ DROP | !("]" | "\"") ~ ANY)* ~ "]" }
WHITESPACE     = _{ " " | "\t" | NEWLINE }
"##]
struct NewickParser;

fn build_graph(input: &str, options: &NewickReadOptions) -> Result<NewickGraph, Report> {
  let tokens = NewickParser::parse(Rule::tree, input).map_err(|error| Report::new(rename_rules(error)))?;
  let mut builder = Builder::new(options);
  for token in tokens {
    builder.consume(token)?;
  }
  builder.finish()
}

fn rename_rules(error: PestError<Rule>) -> PestError<Rule> {
  error.renamed_rules(|rule| {
    match rule {
      Rule::open => "'('",
      Rule::close => "')'",
      Rule::comma => "','",
      Rule::end => "';'",
      Rule::length => "':'",
      Rule::label | Rule::unquoted_label => "label",
      Rule::quoted_label => "quoted label",
      Rule::hybrid_tag => "hybrid tag",
      Rule::number => "branch length",
      Rule::comment | Rule::comment_token | Rule::strict_comment => "comment",
      Rule::rooting | Rule::rooting_value => "rooting comment",
      Rule::EOI => "end of input",
      Rule::tree | Rule::item | Rule::comment_char | Rule::WHITESPACE => "input",
    }
    .to_owned()
  })
}

struct Builder<'i> {
  options: &'i NewickReadOptions,
  graph: NewickGraph,
  hybrids: BTreeMap<(Option<String>, u32), HybridEntry>,
  open: Vec<OpenNode<'i>>,
  current: Slot<'i>,
  seen_token: bool,
  ended: bool,
}

impl<'i> Builder<'i> {
  fn new(options: &'i NewickReadOptions) -> Self {
    Self {
      options,
      graph: NewickGraph::new(),
      hybrids: BTreeMap::new(),
      open: Vec::new(),
      current: Slot::default(),
      seen_token: false,
      ended: false,
    }
  }

  fn consume(&mut self, token: Pair<'i, Rule>) -> Result<(), Report> {
    if self.ended {
      return Ok(());
    }
    match token.as_rule() {
      Rule::rooting => {
        let value = token.into_inner().as_str();
        self.graph.rooted = Some(value.eq_ignore_ascii_case("r"));
        return Ok(());
      },
      Rule::open => {
        if self.current.has_content() {
          return Err(at(
            &token,
            "unexpected '(': a node's '(' must come before its label and branch length",
          ));
        }
        let comments = std::mem::take(&mut self.current.node_comments);
        self.open.push(OpenNode {
          token,
          children: Vec::new(),
          comments,
        });
      },
      Rule::comma => {
        if self.open.is_empty() {
          return Err(at(&token, "unexpected ',' outside of parentheses"));
        }
        let child = self.finish_slot()?;
        if let Some(parent) = self.open.last_mut() {
          parent.children.push(child);
        }
      },
      Rule::close => {
        if self.open.is_empty() {
          return Err(at(&token, "unexpected ')' without a matching '('"));
        }
        let child = self.finish_slot()?;
        if let Some(mut parent) = self.open.pop() {
          parent.children.push(child);
          self.current.children = Some(parent.children);
          self.current.node_comments = parent.comments;
        }
      },
      Rule::label => self.set_label(token)?,
      Rule::length => self.set_length(token)?,
      Rule::comment => self.current.add_comment(token.as_str()),
      Rule::end => self.ended = true,
      _ => return Ok(()),
    }
    self.seen_token = true;
    Ok(())
  }

  fn set_label(&mut self, token: Pair<'i, Rule>) -> Result<(), Report> {
    if self.current.colon_seen {
      return Err(at(
        &token,
        &format!("unexpected label {:?} after the branch length", token.as_str()),
      ));
    }
    if let Some(previous) = &self.current.label {
      return Err(at(
        &token,
        &format!(
          "unexpected label {:?}: the node already has the label {:?}. Quote a label that contains spaces or punctuation",
          token.as_str(),
          previous.text
        ),
      ));
    }
    let text = token.as_str();
    let mut parts = token.into_inner();
    let label = match parts.next() {
      Some(part) if part.as_rule() == Rule::quoted_label => LabelToken {
        text,
        value: unescape_quoted(part.as_str()),
        quoted: true,
        hybrid_tag: parts.next().map(|tag| tag.as_str()),
      },
      _ => LabelToken {
        text,
        value: text.to_owned(),
        quoted: false,
        hybrid_tag: None,
      },
    };
    self.current.label = Some(label);
    Ok(())
  }

  fn set_length(&mut self, token: Pair<'i, Rule>) -> Result<(), Report> {
    if self.current.colon_seen {
      return Err(at(&token, "unexpected ':': the node already has a branch length"));
    }
    self.current.colon_seen = true;
    for part in token.into_inner() {
      match part.as_rule() {
        Rule::number => {
          let number = part
            .as_str()
            .parse::<f64>()
            .wrap_err_with(|| format!("When reading the branch length {:?}", part.as_str()))?;
          self.current.length = Some(number);
        },
        _ => self.current.add_comment(part.as_str()),
      }
    }
    Ok(())
  }

  fn finish(mut self) -> Result<NewickGraph, Report> {
    if let Some(unclosed) = self.open.last() {
      return Err(at(&unclosed.token, "the '(' is never closed"));
    }
    if !self.seen_token {
      return Err(eyre!("The input contains no tree"));
    }
    let (root, root_edge) = self.finish_slot()?;
    let root_node = &mut self.graph.nodes[root];
    root_node.node_attrs.extend(root_edge.branch_attrs);
    root_node.raw_comments.extend(root_edge.raw_comments);
    self.graph.root = root;
    self.graph.validate()?;
    Ok(self.graph)
  }

  fn finish_slot(&mut self) -> Result<(usize, NewickEdgeData), Report> {
    let slot = std::mem::take(&mut self.current);
    let is_internal = slot.children.is_some();
    let parsed = slot
      .label
      .as_ref()
      .map(|label| parse_label(label, is_internal, self.options))
      .transpose()?
      .unwrap_or_default();

    let mut node = NewickNodeData::new();
    node.label = parsed.label;
    node.hybrid = parsed.hybrid.clone();
    for comment in &slot.node_comments {
      classify_comment(comment, &mut node.node_attrs, &mut node.raw_comments);
    }
    let mut edge = NewickEdgeData::new();
    edge.branch_length = slot.length;
    edge.is_acceptor = parsed.is_acceptor;
    for comment in &slot.edge_comments {
      classify_comment(comment, &mut edge.branch_attrs, &mut edge.raw_comments);
    }

    let children = slot.children.unwrap_or_default();
    let idx = match parsed.hybrid {
      Some(hybrid) => self.merge_hybrid(hybrid, node, !children.is_empty())?,
      None => self.graph.add_node(node),
    };
    for (child, child_edge) in children {
      self.graph.add_edge(idx, child, child_edge);
    }
    Ok((idx, edge))
  }

  fn merge_hybrid(&mut self, hybrid: NewickHybrid, node: NewickNodeData, has_children: bool) -> Result<usize, Report> {
    let tag = format!("#{}{}", hybrid.kind.as_deref().unwrap_or(""), hybrid.index);
    let key = (hybrid.kind, hybrid.index);
    let Some(entry) = self.hybrids.get_mut(&key) else {
      let idx = self.graph.add_node(node);
      self.hybrids.insert(key, HybridEntry { idx, has_children });
      return Ok(idx);
    };
    let idx = entry.idx;
    if has_children && entry.has_children {
      return Err(eyre!(
        "The hybrid node {tag} has children in more than one of its occurrences"
      ));
    }
    entry.has_children |= has_children;
    let existing = &mut self.graph.nodes[idx];
    match (&existing.label, node.label) {
      (_, None) => {},
      (None, label) => existing.label = label,
      (Some(previous), Some(label)) if *previous == label => {},
      (Some(_), Some(_)) => {
        return Err(eyre!(
          "The occurrences of the hybrid node {tag} have different labels, at {}",
          describe_node(&self.graph, idx)
        ));
      },
    }
    let existing = &mut self.graph.nodes[idx];
    existing.node_attrs.extend(node.node_attrs);
    existing.raw_comments.extend(node.raw_comments);
    Ok(idx)
  }
}

#[derive(Default)]
struct Slot<'i> {
  label: Option<LabelToken<'i>>,
  colon_seen: bool,
  length: Option<f64>,
  node_comments: Vec<&'i str>,
  edge_comments: Vec<&'i str>,
  children: Option<Vec<(usize, NewickEdgeData)>>,
}

impl<'i> Slot<'i> {
  fn has_content(&self) -> bool {
    self.label.is_some() || self.colon_seen || self.children.is_some()
  }

  fn add_comment(&mut self, comment: &'i str) {
    if self.colon_seen {
      self.edge_comments.push(comment);
    } else {
      self.node_comments.push(comment);
    }
  }
}

struct OpenNode<'i> {
  token: Pair<'i, Rule>,
  children: Vec<(usize, NewickEdgeData)>,
  comments: Vec<&'i str>,
}

struct LabelToken<'i> {
  text: &'i str,
  value: String,
  quoted: bool,
  hybrid_tag: Option<&'i str>,
}

struct HybridEntry {
  idx: usize,
  has_children: bool,
}

#[derive(Default)]
struct ParsedLabel {
  label: Option<NewickLabel>,
  hybrid: Option<NewickHybrid>,
  is_acceptor: bool,
}

fn parse_label(token: &LabelToken<'_>, is_internal: bool, options: &NewickReadOptions) -> Result<ParsedLabel, Report> {
  let (name, hybrid_tag) = match (token.quoted, token.hybrid_tag) {
    (true, Some(tag)) if options.enewick => (token.value.clone(), Some(tag)),
    (true, Some(tag)) => (format!("{}{tag}", token.value), None),
    (true, None) => (token.value.clone(), None),
    (false, _) if options.enewick => split_unquoted_hybrid(&token.value),
    (false, _) => (token.value.clone(), None),
  };
  let (hybrid, is_acceptor) = match hybrid_tag {
    Some(tag) => {
      let (hybrid, is_acceptor) = parse_hybrid_tag(tag, token.text)?;
      (Some(hybrid), is_acceptor)
    },
    None => (None, false),
  };
  let label = if name.is_empty() && hybrid.is_some() {
    None
  } else if is_internal {
    Some(parse_support(&name).map_or(NewickLabel::Name(name), NewickLabel::Support))
  } else {
    Some(NewickLabel::Name(name))
  };
  Ok(ParsedLabel {
    label,
    hybrid,
    is_acceptor,
  })
}

#[expect(
  clippy::string_slice,
  reason = "the indices come from regex matches on the same string"
)]
fn split_unquoted_hybrid(text: &str) -> (String, Option<&str>) {
  match regex!(r"^(.*?)(##?[A-Za-z]*[0-9]+)$").captures(text) {
    Some(captures) => match (captures.get(1), captures.get(2)) {
      (Some(name), Some(tag)) => (name.as_str().to_owned(), Some(&text[tag.start()..tag.end()])),
      _ => (text.to_owned(), None),
    },
    None => (text.to_owned(), None),
  }
}

fn parse_hybrid_tag(tag: &str, label: &str) -> Result<(NewickHybrid, bool), Report> {
  let is_acceptor = tag.starts_with("##");
  let body = tag.trim_start_matches('#');
  let digits_start = body.find(|c: char| c.is_ascii_digit()).unwrap_or(body.len());
  let (kind, index) = body.split_at(digits_start);
  let index = index
    .parse::<u32>()
    .wrap_err_with(|| format!("When parsing the hybrid node index in '{label}'"))?;
  let kind = (!kind.is_empty()).then(|| kind.to_owned());
  Ok((NewickHybrid { kind, index }, is_acceptor))
}

#[cfg_attr(
  dylint_lib = "treetime_lints",
  expect(
    result_defaulted,
    reason = "an internal label that is not a number is a node name, not a malformed support value"
  )
)]
fn parse_support(label: &str) -> Option<f64> {
  label.parse::<f64>().ok()
}

fn unescape_quoted(quoted: &str) -> String {
  quoted
    .strip_prefix('\'')
    .and_then(|inner| inner.strip_suffix('\''))
    .unwrap_or(quoted)
    .replace("''", "'")
}

fn at(token: &Pair<'_, Rule>, message: &str) -> Report {
  let (line, column) = token.line_col();
  eyre!("At line {line}, column {column}: {message}")
}
