//! Newick parser with automatic annotation dialect detection.

use crate::types::{NewickEdgeData, NewickGraph, NewickHybrid, NewickNodeData, NewickValue};
use eyre::{Report, WrapErr};
use pest::Parser;
use pest_derive::Parser;
use regex::Regex;
use std::collections::{BTreeMap, HashMap};
use std::io::Read;
use std::sync::LazyLock;

#[derive(Parser)]
#[grammar_inline = r#"
tree        = _{ SOI ~ rooting? ~ root_branch ~ ";" ~ EOI }
root_branch =  { subtree ~ comment* ~ (":" ~ comment* ~ number)? ~ comment* }
subtree     =  { leaf | internal }
internal =  { "(" ~ branch_list ~ ")" ~ label? ~ comment* }
leaf     =  { label ~ comment* }

branch_list = { branch ~ ("," ~ branch)* }
branch      = { subtree? ~ comment* ~ (":" ~ comment* ~ number)? ~ comment* }

rooting        = { "[&" ~ rooting_value ~ "]" }
rooting_value  = @{ ^"r" | ^"u" }

label          = { quoted_label | unquoted_label }
quoted_label   = @{ "'" ~ (!"'" ~ ANY | "''")* ~ "'" }
unquoted_label = @{ safe_char+ }
safe_char      = _{ !("(" | ")" | "[" | "]" | "," | ";" | ":" | WHITESPACE) ~ ANY }
number  = @{ "-"? ~ ASCII_DIGIT+ ~ ("." ~ ASCII_DIGIT*)? ~ (^"e" ~ ("+" | "-")? ~ ASCII_DIGIT+)? }

comment       = @{ "[" ~ comment_inner ~ "]" }
comment_inner = _{ (!"[" ~ !"]" ~ ANY | "[" ~ comment_inner ~ "]")* }

WHITESPACE = _{ " " | "\t" | NEWLINE }
"#]
struct NewickParser;

static HYBRID_RE: LazyLock<Regex> =
  LazyLock::new(|| Regex::new(r"^(.*?)(##|#)([A-Za-z]*)(\d+)$").expect("hybrid regex"));

/// Parse a Newick string into a [`NewickGraph`], auto-detecting BEAST, NHX, and plain comment dialects.
pub fn newick_from_string(input: &str) -> Result<NewickGraph, Report> {
  let pairs = NewickParser::parse(Rule::tree, input).wrap_err("Failed to parse Newick string")?;

  let mut graph = NewickGraph::new();
  let mut hybrid_map: HashMap<(Option<String>, u32), usize> = HashMap::new();

  let mut rooted = None;
  let mut root = None;

  for pair in pairs {
    match pair.as_rule() {
      Rule::rooting => {
        let inner = pair.into_inner().next().expect("rooting_value");
        rooted = Some(inner.as_str().eq_ignore_ascii_case("r"));
      },
      Rule::root_branch => {
        root = Some(visit_root_branch(pair, &mut graph, &mut hybrid_map)?);
      },
      _ => {},
    }
  }

  graph.root = root.expect("tree must have a root subtree");
  graph.rooted = rooted;

  Ok(graph)
}

/// Parse Newick from an `impl Read` source.
pub fn newick_from_reader(mut reader: impl Read) -> Result<NewickGraph, Report> {
  let mut input = String::new();
  reader.read_to_string(&mut input)?;
  newick_from_string(&input)
}

fn visit_root_branch(
  pair: pest::iterators::Pair<Rule>,
  graph: &mut NewickGraph,
  hybrid_map: &mut HashMap<(Option<String>, u32), usize>,
) -> Result<usize, Report> {
  let mut root_idx = None;
  let mut root_attrs = BTreeMap::new();
  let mut root_raw = Vec::new();

  for inner in pair.into_inner() {
    match inner.as_rule() {
      Rule::subtree => {
        let (idx, _acceptor) = visit_subtree(inner, graph, hybrid_map)?;
        root_idx = Some(idx);
      },
      Rule::number => {
        // Root branch length discarded (no parent edge to attach it to)
      },
      Rule::comment => {
        classify_comment(inner.as_str(), &mut root_attrs, &mut root_raw);
      },
      _ => {},
    }
  }

  let idx = root_idx.expect("root_branch must contain a subtree");
  let root_node = &mut graph.nodes[idx];
  root_node.node_attrs.extend(root_attrs);
  root_node.raw_comments.extend(root_raw);
  Ok(idx)
}

fn visit_subtree(
  pair: pest::iterators::Pair<Rule>,
  graph: &mut NewickGraph,
  hybrid_map: &mut HashMap<(Option<String>, u32), usize>,
) -> Result<(usize, bool), Report> {
  let inner = pair.into_inner().next().expect("subtree must have leaf or internal");
  match inner.as_rule() {
    Rule::leaf => visit_leaf(inner, graph, hybrid_map),
    Rule::internal => visit_internal(inner, graph, hybrid_map),
    _ => unreachable!(),
  }
}

fn visit_leaf(
  pair: pest::iterators::Pair<Rule>,
  graph: &mut NewickGraph,
  hybrid_map: &mut HashMap<(Option<String>, u32), usize>,
) -> Result<(usize, bool), Report> {
  let mut name = None;
  let mut raw_comments = Vec::new();
  let mut node_attrs = BTreeMap::new();

  for inner in pair.into_inner() {
    match inner.as_rule() {
      Rule::label => {
        name = Some(parse_label(inner));
      },
      Rule::comment => {
        classify_comment(inner.as_str(), &mut node_attrs, &mut raw_comments);
      },
      _ => {},
    }
  }

  let extraction = extract_hybrid(name);
  let node_data = NewickNodeData {
    name: extraction.clean_name,
    confidence: None,
    node_attrs,
    raw_comments,
    hybrid: extraction.hybrid.clone(),
    children: Vec::new(),
  };

  let idx = resolve_hybrid_node(node_data, extraction.hybrid.as_ref(), graph, hybrid_map);
  Ok((idx, extraction.is_acceptor))
}

fn visit_internal(
  pair: pest::iterators::Pair<Rule>,
  graph: &mut NewickGraph,
  hybrid_map: &mut HashMap<(Option<String>, u32), usize>,
) -> Result<(usize, bool), Report> {
  let mut name = None;
  let mut raw_comments = Vec::new();
  let mut node_attrs = BTreeMap::new();
  let mut branch_children = Vec::new();

  for inner in pair.into_inner() {
    match inner.as_rule() {
      Rule::branch_list => {
        for branch in inner.into_inner() {
          if branch.as_rule() == Rule::branch {
            branch_children.push(visit_branch(branch, graph, hybrid_map)?);
          }
        }
      },
      Rule::label => {
        name = Some(parse_label(inner));
      },
      Rule::comment => {
        classify_comment(inner.as_str(), &mut node_attrs, &mut raw_comments);
      },
      _ => {},
    }
  }

  let extraction = extract_hybrid(name);

  // Biopython heuristic: bare numeric label on an internal node is branch support,
  // not a taxon name. Try parse as float; on success, store as confidence and clear name.
  let (final_name, confidence) = match extraction.clean_name {
    Some(ref label) if label.parse::<f64>().is_ok() => (None, label.parse::<f64>().ok()),
    other => (other, None),
  };

  let node_data = NewickNodeData {
    name: final_name,
    confidence,
    node_attrs,
    raw_comments,
    hybrid: extraction.hybrid.clone(),
    children: Vec::new(),
  };

  let node_idx = resolve_hybrid_node(node_data, extraction.hybrid.as_ref(), graph, hybrid_map);

  for (child_idx, edge_data) in branch_children {
    graph.add_edge(node_idx, child_idx, edge_data);
  }

  Ok((node_idx, extraction.is_acceptor))
}

fn visit_branch(
  pair: pest::iterators::Pair<Rule>,
  graph: &mut NewickGraph,
  hybrid_map: &mut HashMap<(Option<String>, u32), usize>,
) -> Result<(usize, NewickEdgeData), Report> {
  let mut child_idx = None;
  let mut branch_length = None;
  let mut seen_number = false;

  let mut is_acceptor = false;
  let mut raw_comments: Vec<String> = Vec::new();
  let mut pre_colon_attrs = BTreeMap::new();
  let mut post_colon_attrs = BTreeMap::new();

  for inner in pair.into_inner() {
    match inner.as_rule() {
      Rule::subtree => {
        let (idx, acceptor) = visit_subtree(inner, graph, hybrid_map)?;
        child_idx = Some(idx);
        is_acceptor = acceptor;
      },
      Rule::number => {
        branch_length = Some(
          inner
            .as_str()
            .parse::<f64>()
            .wrap_err_with(|| format!("Invalid branch length: '{}'", inner.as_str()))?,
        );
        seen_number = true;
      },
      Rule::comment => {
        if seen_number {
          classify_comment(inner.as_str(), &mut post_colon_attrs, &mut raw_comments);
        } else {
          classify_comment(inner.as_str(), &mut pre_colon_attrs, &mut raw_comments);
        }
      },
      _ => {},
    }
  }

  // If no subtree was present, create an anonymous leaf node
  let child_idx = child_idx.unwrap_or_else(|| graph.add_node(NewickNodeData::new()));

  let mut branch_attrs = pre_colon_attrs;
  branch_attrs.extend(post_colon_attrs);

  let edge_data = NewickEdgeData {
    branch_length,
    branch_attrs,
    raw_comments,
    is_acceptor,
  };

  Ok((child_idx, edge_data))
}

fn parse_label(pair: pest::iterators::Pair<Rule>) -> String {
  let inner = pair.into_inner().next().expect("label must have content");
  match inner.as_rule() {
    Rule::quoted_label => {
      let raw = inner.as_str();
      raw
        .strip_prefix('\'')
        .and_then(|s| s.strip_suffix('\''))
        .unwrap_or(raw)
        .replace("''", "'")
    },
    Rule::unquoted_label => inner.as_str().to_owned(),
    _ => unreachable!(),
  }
}

struct HybridExtraction {
  clean_name: Option<String>,
  hybrid: Option<NewickHybrid>,
  is_acceptor: bool,
}

fn extract_hybrid(name: Option<String>) -> HybridExtraction {
  let Some(name) = name else {
    return HybridExtraction {
      clean_name: None,
      hybrid: None,
      is_acceptor: false,
    };
  };

  let Some(caps) = HYBRID_RE.captures(&name) else {
    return HybridExtraction {
      clean_name: Some(name),
      hybrid: None,
      is_acceptor: false,
    };
  };

  let prefix = caps.get(1).map_or("", |m| m.as_str());
  let hash_marker = caps.get(2).map_or("", |m| m.as_str());
  let kind_str = caps.get(3).map_or("", |m| m.as_str());
  let index_str = caps.get(4).map_or("0", |m| m.as_str());

  let clean_name = if prefix.is_empty() {
    None
  } else {
    Some(prefix.to_owned())
  };
  let kind = if kind_str.is_empty() {
    None
  } else {
    Some(kind_str.to_owned())
  };
  let index = index_str.parse::<u32>().unwrap_or(0);
  let is_acceptor = hash_marker == "##";

  HybridExtraction {
    clean_name,
    hybrid: Some(NewickHybrid { kind, index }),
    is_acceptor,
  }
}

fn resolve_hybrid_node(
  node_data: NewickNodeData,
  hybrid: Option<&NewickHybrid>,
  graph: &mut NewickGraph,
  hybrid_map: &mut HashMap<(Option<String>, u32), usize>,
) -> usize {
  let Some(h) = hybrid else {
    return graph.add_node(node_data);
  };

  let key = (h.kind.clone(), h.index);
  if let Some(&existing_idx) = hybrid_map.get(&key) {
    let existing = &mut graph.nodes[existing_idx];
    existing.node_attrs.extend(node_data.node_attrs);
    existing.raw_comments.extend(node_data.raw_comments);
    existing_idx
  } else {
    let idx = graph.add_node(node_data);
    hybrid_map.insert(key, idx);
    idx
  }
}

fn classify_comment(comment_text: &str, attrs: &mut BTreeMap<String, NewickValue>, raw: &mut Vec<String>) {
  let inner = comment_text
    .strip_prefix('[')
    .and_then(|s| s.strip_suffix(']'))
    .unwrap_or(comment_text);

  if let Some(nhx_body) = inner.strip_prefix("&&NHX") {
    parse_nhx_attrs(nhx_body, attrs);
  } else if let Some(beast_body) = inner.strip_prefix('&') {
    parse_beast_attrs(beast_body, attrs);
  } else {
    raw.push(comment_text.to_owned());
  }
}

fn parse_nhx_attrs(body: &str, attrs: &mut BTreeMap<String, NewickValue>) {
  // NHX format: :key=value:key=value...
  // Leading colon before first tag
  for tag in body.split(':') {
    let tag = tag.trim();
    if tag.is_empty() {
      continue;
    }
    if let Some((key, value)) = tag.split_once('=') {
      attrs.insert(key.to_owned(), NewickValue::String(value.to_owned()));
    }
  }
}

fn parse_beast_attrs(body: &str, attrs: &mut BTreeMap<String, NewickValue>) {
  // BEAST format: key=value,key=value,...
  // Need to handle {array} values and "quoted strings"
  let pairs = split_beast_pairs(body);
  for pair_str in pairs {
    let pair_str = pair_str.trim();
    if pair_str.is_empty() {
      continue;
    }
    if let Some((key, value_str)) = pair_str.split_once('=') {
      let key = strip_beast_quotes(key.trim());
      let value_str = value_str.trim();
      attrs.insert(key, parse_beast_value(value_str));
    } else {
      // Bare key = boolean TRUE
      attrs.insert(strip_beast_quotes(pair_str), NewickValue::Boolean(true));
    }
  }
}

fn split_beast_pairs(body: &str) -> Vec<String> {
  let mut pairs = Vec::new();
  let mut current = String::new();
  let mut brace_depth = 0_u32;
  let mut in_quote = false;

  for ch in body.chars() {
    match ch {
      '"' if brace_depth == 0 => {
        in_quote = !in_quote;
        current.push(ch);
      },
      '{' if !in_quote => {
        brace_depth += 1;
        current.push(ch);
      },
      '}' if !in_quote => {
        brace_depth = brace_depth.saturating_sub(1);
        current.push(ch);
      },
      ',' if brace_depth == 0 && !in_quote => {
        pairs.push(std::mem::take(&mut current));
      },
      _ => {
        current.push(ch);
      },
    }
  }
  if !current.is_empty() {
    pairs.push(current);
  }
  pairs
}

fn strip_beast_quotes(s: &str) -> String {
  if let Some(inner) = s.strip_prefix('"').and_then(|i| i.strip_suffix('"')) {
    return inner.replace("\"\"", "\"");
  }
  if let Some(inner) = s.strip_prefix('\'').and_then(|i| i.strip_suffix('\'')) {
    return inner.replace("''", "'");
  }
  s.to_owned()
}

fn parse_beast_value(s: &str) -> NewickValue {
  let s = s.trim();

  if let Some(inner) = s.strip_prefix('{').and_then(|s| s.strip_suffix('}')) {
    let elements: Vec<NewickValue> = split_beast_pairs(inner)
      .iter()
      .map(|e| parse_beast_value(e.trim()))
      .collect();
    return NewickValue::Array(elements);
  }

  // Boolean: TRUE/FALSE (case-insensitive)
  if s.eq_ignore_ascii_case("true") {
    return NewickValue::Boolean(true);
  }
  if s.eq_ignore_ascii_case("false") {
    return NewickValue::Boolean(false);
  }

  // Number: starts with digit or '-', must be finite
  if s.starts_with(|c: char| c.is_ascii_digit() || c == '-') {
    if let Ok(n) = s.parse::<f64>() {
      if n.is_finite() {
        return NewickValue::Number(n);
      }
    }
  }

  // String (strip quotes if present)
  NewickValue::String(strip_beast_quotes(s))
}
