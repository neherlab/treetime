//! Newick writer with configurable annotation style.

use crate::types::{NewickEdgeData, NewickGraph, NewickNodeData, NewickValue, NewickWriteOptions, NwkStyle};
use eyre::Report;
use std::collections::{BTreeMap, HashSet};
use std::io::Write;

/// Write a [`NewickGraph`] as a Newick string.
pub fn newick_to_string(graph: &NewickGraph, options: &NewickWriteOptions) -> Result<String, Report> {
  let mut buf = Vec::new();
  newick_to_writer(&mut buf, graph, options)?;
  Ok(String::from_utf8(buf).expect("Newick output should be valid UTF-8"))
}

/// Write a [`NewickGraph`] to an `impl Write` destination.
pub fn newick_to_writer(
  writer: &mut impl Write,
  graph: &NewickGraph,
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  if let Some(rooted) = graph.rooted {
    if rooted {
      write!(writer, "[&R]")?;
    } else {
      write!(writer, "[&U]")?;
    }
  }

  let mut visited_hybrids = HashSet::new();
  write_subtree(writer, graph, graph.root, false, options, &mut visited_hybrids)?;
  write!(writer, ";")?;

  Ok(())
}

fn write_subtree(
  writer: &mut impl Write,
  graph: &NewickGraph,
  node_idx: usize,
  is_acceptor: bool,
  options: &NewickWriteOptions,
  visited_hybrids: &mut HashSet<usize>,
) -> Result<(), Report> {
  let node = &graph.nodes[node_idx];

  // For hybrid nodes: if already visited, write marker only (no subtree)
  if node.hybrid.is_some() && !visited_hybrids.insert(node_idx) {
    write_name(writer, node, is_acceptor)?;
    return Ok(());
  }

  let children = &node.children;

  if !children.is_empty() {
    write!(writer, "(")?;
    for (i, &edge_idx) in children.iter().enumerate() {
      if i > 0 {
        write!(writer, ",")?;
      }
      let edge = &graph.edges[edge_idx];
      write_subtree(
        writer,
        graph,
        edge.child,
        edge.data.is_acceptor,
        options,
        visited_hybrids,
      )?;
      write_edge(writer, &edge.data, options)?;
    }
    write!(writer, ")")?;
  }

  write_name(writer, node, is_acceptor)?;
  write_node_annotations(writer, node, options)?;

  Ok(())
}

fn write_name(writer: &mut impl Write, node: &NewickNodeData, is_acceptor: bool) -> Result<(), Report> {
  if node.hybrid.is_none() {
    if let Some(name) = &node.name {
      return write_label(writer, name);
    }
    return Ok(());
  }

  let mut label = String::new();
  if let Some(name) = &node.name {
    label.push_str(name);
  }
  if let Some(hybrid) = &node.hybrid {
    label.push('#');
    if is_acceptor {
      label.push('#');
    }
    if let Some(kind) = &hybrid.kind {
      label.push_str(kind);
    }
    label.push_str(&hybrid.index.to_string());
  }
  if !label.is_empty() {
    write_label(writer, &label)?;
  }
  Ok(())
}

fn write_label(writer: &mut impl Write, label: &str) -> Result<(), Report> {
  if needs_quoting(label) {
    write!(writer, "'")?;
    for part in label.split_inclusive('\'') {
      if let Some(prefix) = part.strip_suffix('\'') {
        write!(writer, "{prefix}''")?;
      } else {
        write!(writer, "{part}")?;
      }
    }
    write!(writer, "'")?;
  } else {
    write!(writer, "{label}")?;
  }
  Ok(())
}

fn needs_quoting(name: &str) -> bool {
  name.contains(|c: char| matches!(c, '(' | ')' | '[' | ']' | ',' | ';' | ':' | '\'') || c.is_whitespace())
}

fn write_node_annotations(
  writer: &mut impl Write,
  node: &NewickNodeData,
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  match options.style {
    NwkStyle::Plain => {},
    NwkStyle::Beast => {
      write_beast_attrs(writer, &node.node_attrs)?;
      write_raw_comments(writer, &node.raw_comments)?;
    },
    NwkStyle::Nhx => {
      write_nhx_attrs(writer, &node.node_attrs)?;
      write_raw_comments(writer, &node.raw_comments)?;
    },
  }
  Ok(())
}

fn write_edge(writer: &mut impl Write, edge: &NewickEdgeData, options: &NewickWriteOptions) -> Result<(), Report> {
  let has_length = edge.branch_length.is_some();
  let has_branch_attrs = !edge.branch_attrs.is_empty();
  let has_raw = !edge.raw_comments.is_empty();

  if !has_length && !has_branch_attrs && !has_raw {
    return Ok(());
  }

  if has_length {
    write!(writer, ":")?;

    // BEAST2 canonical: branch attrs between `:` and number
    if options.style == NwkStyle::Beast {
      write_beast_attrs(writer, &edge.branch_attrs)?;
    }

    write!(writer, "{}", format_float(edge.branch_length.unwrap_or(0.0), options))?;

    // NHX and MrBayes: branch attrs after number
    if options.style == NwkStyle::Nhx {
      write_nhx_attrs(writer, &edge.branch_attrs)?;
    }

    // Raw comments after number for both Beast and Nhx
    if options.style != NwkStyle::Plain {
      write_raw_comments(writer, &edge.raw_comments)?;
    }
  } else if options.style != NwkStyle::Plain && (has_branch_attrs || has_raw) {
    // No length: emit attrs and raw comments without `:` prefix.
    // These will be parsed back as pre-colon branch comments.
    match options.style {
      NwkStyle::Beast => write_beast_attrs(writer, &edge.branch_attrs)?,
      NwkStyle::Nhx => write_nhx_attrs(writer, &edge.branch_attrs)?,
      NwkStyle::Plain => {},
    }
    write_raw_comments(writer, &edge.raw_comments)?;
  }

  Ok(())
}

fn write_beast_attrs(writer: &mut impl Write, attrs: &BTreeMap<String, NewickValue>) -> Result<(), Report> {
  if attrs.is_empty() {
    return Ok(());
  }
  write!(writer, "[&")?;
  for (i, (key, value)) in attrs.iter().enumerate() {
    if i > 0 {
      write!(writer, ",")?;
    }
    if needs_beast_quoting(key) {
      let escaped = key.replace('"', "\"\"");
      write!(writer, "\"{escaped}\"=")?;
    } else {
      write!(writer, "{key}=")?;
    }
    write_beast_value(writer, value)?;
  }
  write!(writer, "]")?;
  Ok(())
}

fn write_beast_value(writer: &mut impl Write, value: &NewickValue) -> Result<(), Report> {
  match value {
    NewickValue::Boolean(b) => {
      if *b {
        write!(writer, "TRUE")?;
      } else {
        write!(writer, "FALSE")?;
      }
    },
    NewickValue::Number(n) => {
      write!(writer, "{n}")?;
    },
    NewickValue::String(s) => {
      if needs_beast_quoting(s) {
        let escaped = s.replace('"', "\"\"");
        write!(writer, "\"{escaped}\"")?;
      } else {
        write!(writer, "{s}")?;
      }
    },
    NewickValue::Array(arr) => {
      write!(writer, "{{")?;
      for (i, elem) in arr.iter().enumerate() {
        if i > 0 {
          write!(writer, ",")?;
        }
        write_beast_value(writer, elem)?;
      }
      write!(writer, "}}")?;
    },
  }
  Ok(())
}

fn needs_beast_quoting(s: &str) -> bool {
  if s.contains(|c: char| matches!(c, ',' | '=' | ']' | '[' | '\'' | '"') || c.is_whitespace()) {
    return true;
  }
  if s.eq_ignore_ascii_case("true") || s.eq_ignore_ascii_case("false") {
    return true;
  }
  if s.starts_with(|c: char| c.is_ascii_digit() || c == '-') && s.parse::<f64>().is_ok() {
    return true;
  }
  false
}

fn write_nhx_attrs(writer: &mut impl Write, attrs: &BTreeMap<String, NewickValue>) -> Result<(), Report> {
  if attrs.is_empty() {
    return Ok(());
  }
  write!(writer, "[&&NHX")?;
  for (key, value) in attrs {
    write!(writer, ":{key}=")?;
    write_nhx_value(writer, value)?;
  }
  write!(writer, "]")?;
  Ok(())
}

fn write_nhx_value(writer: &mut impl Write, value: &NewickValue) -> Result<(), Report> {
  match value {
    NewickValue::Boolean(b) => write!(writer, "{b}")?,
    NewickValue::Number(n) => write!(writer, "{n}")?,
    NewickValue::String(s) => {
      if s.contains(|c: char| matches!(c, ':' | '=' | ']')) {
        return Err(eyre::eyre!(
          "NHX cannot represent value containing reserved character (':', '=', or ']'): {s}"
        ));
      }
      write!(writer, "{s}")?;
    },
    NewickValue::Array(arr) => {
      for (i, elem) in arr.iter().enumerate() {
        if i > 0 {
          write!(writer, ">")?;
        }
        write_nhx_value(writer, elem)?;
      }
    },
  }
  Ok(())
}

fn write_raw_comments(writer: &mut impl Write, comments: &[String]) -> Result<(), Report> {
  for comment in comments {
    write!(writer, "{comment}")?;
  }
  Ok(())
}

fn format_float(value: f64, options: &NewickWriteOptions) -> String {
  let mut config = pretty_dtoa::FmtFloatConfig::default()
    .add_point_zero(false)
    .radix_point('.');

  if let Some(sig) = options.significant_digits {
    config = config.max_significant_digits(sig);
  }
  if let Some(dec) = options.decimal_digits {
    config = config.max_decimal_digits(dec);
  }

  pretty_dtoa::dtoa(value, config)
}
