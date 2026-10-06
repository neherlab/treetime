use crate::annotation::{write_beast_attrs, write_nhx_attrs, write_raw_comments};
use crate::number::{format_number, format_shortest};
use crate::types::{NewickEdgeData, NewickGraph, NewickLabel, NewickNodeData, NewickWriteOptions, NwkStyle};
use crate::validate::describe_node;
use eyre::{Report, WrapErr, eyre};
use std::io;

pub fn newick_to_writer(
  writer: &mut impl io::Write,
  graph: &NewickGraph,
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  write_newick(writer, graph, options).wrap_err("When writing Newick")
}

pub fn newick_to_string(graph: &NewickGraph, options: &NewickWriteOptions) -> Result<String, Report> {
  let mut buffer = Vec::new();
  newick_to_writer(&mut buffer, graph, options)?;
  Ok(String::from_utf8(buffer)?)
}

pub fn write_label(writer: &mut impl io::Write, label: &str) -> Result<(), Report> {
  if needs_quoting(label) {
    write!(writer, "'{}'", label.replace('\'', "''"))?;
  } else {
    writer.write_all(label.as_bytes())?;
  }
  Ok(())
}

pub fn needs_quoting(name: &str) -> bool {
  name.is_empty()
    || name.contains(|c: char| matches!(c, '(' | ')' | '[' | ']' | ',' | ';' | ':' | '\'' | '#') || c.is_whitespace())
}

pub(crate) fn write_newick(
  writer: &mut impl io::Write,
  graph: &NewickGraph,
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  graph.validate()?;
  if options.significant_digits == Some(0) {
    return Err(eyre!("The number of significant digits must be at least 1"));
  }
  match graph.rooted {
    Some(true) => writer.write_all(b"[&R]")?,
    Some(false) => writer.write_all(b"[&U]")?,
    None => {},
  }

  let mut defined = vec![false; graph.nodes.len()];
  defined[graph.root] = true;
  let mut stack = vec![Frame {
    node: graph.root,
    edge: None,
    next_child: 0,
    is_definition: true,
  }];
  while let Some(frame) = stack.last_mut() {
    let node = &graph.nodes[frame.node];
    let children: &[usize] = if frame.is_definition { &node.children } else { &[] };
    if let Some(&edge_idx) = children.get(frame.next_child) {
      writer.write_all(if frame.next_child == 0 { b"(" } else { b"," })?;
      frame.next_child += 1;
      let child = graph.edges[edge_idx].child;
      let is_definition = graph.nodes[child].hybrid.is_none() || !std::mem::replace(&mut defined[child], true);
      stack.push(Frame {
        node: child,
        edge: Some(edge_idx),
        next_child: 0,
        is_definition,
      });
      continue;
    }
    let Frame {
      node: node_idx,
      edge,
      next_child,
      is_definition,
    } = *frame;
    stack.pop();
    if next_child > 0 {
      writer.write_all(b")")?;
    }
    let is_acceptor = edge.is_some_and(|edge_idx| graph.edges[edge_idx].data.is_acceptor);
    write_node_label(writer, graph, node_idx, is_definition, is_acceptor)
      .wrap_err_with(|| format!("When writing the label of {}", describe_node(graph, node_idx)))?;
    if is_definition {
      write_node_annotations(writer, node, options)?;
    }
    if let Some(edge_idx) = edge {
      write_edge(writer, &graph.edges[edge_idx].data, options)
        .wrap_err_with(|| format!("When writing the branch above {}", describe_node(graph, node_idx)))?;
    }
  }
  writer.write_all(b";")?;
  Ok(())
}

#[derive(Clone, Copy)]
struct Frame {
  node: usize,
  edge: Option<usize>,
  next_child: usize,
  is_definition: bool,
}

fn write_node_label(
  writer: &mut impl io::Write,
  graph: &NewickGraph,
  node_idx: usize,
  is_definition: bool,
  is_acceptor: bool,
) -> Result<(), Report> {
  let node = &graph.nodes[node_idx];
  let is_internal = !node.children.is_empty();
  match &node.label {
    Some(NewickLabel::Name(name)) => {
      if is_internal && name.parse::<f64>().is_ok() {
        return Err(eyre!(
          "The internal node name {name:?} would be read back as a support value"
        ));
      }
      write_label(writer, name)?;
    },
    Some(NewickLabel::Support(support)) if is_definition => {
      if !is_internal {
        return Err(eyre!(
          "A leaf cannot carry a support value, because Newick reads a leaf label as a name"
        ));
      }
      let text = if support.is_finite() {
        format_shortest(*support)?
      } else {
        support.to_string()
      };
      writer.write_all(text.as_bytes())?;
    },
    Some(NewickLabel::Support(_)) | None => {},
  }
  if let Some(hybrid) = &node.hybrid {
    let marker = if is_acceptor { "##" } else { "#" };
    write!(
      writer,
      "{marker}{}{}",
      hybrid.kind.as_deref().unwrap_or(""),
      hybrid.index
    )?;
  }
  Ok(())
}

fn write_node_annotations(
  writer: &mut impl io::Write,
  node: &NewickNodeData,
  options: &NewickWriteOptions,
) -> Result<(), Report> {
  let attrs = node.node_attrs.iter().map(|(key, value)| (key.as_str(), value));
  match options.style {
    NwkStyle::Plain => return Ok(()),
    NwkStyle::Beast => write_beast_attrs(writer, attrs)?,
    NwkStyle::Nhx => write_nhx_attrs(writer, attrs)?,
  }
  write_raw_comments(writer, &node.raw_comments)
}

fn write_edge(writer: &mut impl io::Write, edge: &NewickEdgeData, options: &NewickWriteOptions) -> Result<(), Report> {
  let has_comments =
    options.style != NwkStyle::Plain && (!edge.branch_attrs.is_empty() || !edge.raw_comments.is_empty());
  if edge.branch_length.is_none() && !has_comments {
    return Ok(());
  }
  writer.write_all(b":")?;
  let attrs = edge.branch_attrs.iter().map(|(key, value)| (key.as_str(), value));
  if options.style == NwkStyle::Beast {
    write_beast_attrs(writer, attrs.clone())?;
  }
  if let Some(length) = edge.branch_length {
    let text = format_number(length, options.significant_digits, options.decimal_digits)?;
    writer.write_all(text.as_bytes())?;
  }
  if options.style == NwkStyle::Nhx {
    write_nhx_attrs(writer, attrs)?;
  }
  if options.style != NwkStyle::Plain {
    write_raw_comments(writer, &edge.raw_comments)?;
  }
  Ok(())
}
