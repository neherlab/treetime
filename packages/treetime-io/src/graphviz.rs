use crate::nwk::{NwkWriteOptions, format_weight};
use eyre::{Report, WrapErr};
use itertools::{Itertools, iproduct};
use std::collections::BTreeMap;
use std::fmt::Write;
use std::io;
use std::path::Path;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, Node};
use treetime_utils::io::file::write_file_with;
use treetime_utils::make_internal_report;

const FAKE_EDGE_BASE_WEIGHT: usize = 1000;
const FAKE_EDGE_WEIGHT_STEP: usize = 100;

pub fn graphviz_write_file(
  filepath: impl AsRef<Path>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  write_file_with(filepath, |writer| graphviz_write(writer, graph, names, weights))
}

pub fn graphviz_write(
  mut writer: impl io::Write,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  let mut text = String::new();
  write_graphviz_text(&mut text, graph, names, weights)?;
  writer.write_all(text.as_bytes()).wrap_err("When writing Graphviz")
}

fn write_graphviz_text(
  mut writer: impl Write,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  write!(
    writer,
    r#"digraph Phylogeny {{
  graph [rankdir=LR, overlap=scale, splines=ortho, nodesep=1.0, ordering=out];
  edge  [overlap=scale];
  node  [shape=box];
"#
  )?;
  print_nodes(graph, names, &mut writer)?;
  writeln!(writer)?;
  print_edges(graph, weights, &mut writer)?;
  writeln!(writer, "}}")?;
  Ok(())
}

fn print_nodes<W>(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>, mut writer: W) -> Result<(), Report>
where
  W: Write,
{
  writeln!(writer, "\n  subgraph roots {{")?;
  let roots = graph.get_roots().collect::<Vec<_>>();
  for &node in &roots {
    print_node(&mut writer, node, names)?;
  }
  print_fake_edges(&mut writer, &roots.iter().map(|node| node.key()).collect_vec())?;

  writeln!(writer, "  }}\n\n  subgraph internals {{")?;
  for node in graph.get_internal_nodes().filter(|node| !node.is_root()) {
    print_node(&mut writer, node, names)?;
  }

  writeln!(writer, "  }}\n\n  subgraph leaves {{")?;
  let leaves = graph.get_leaves().filter(|node| !node.is_root()).collect::<Vec<_>>();
  for &node in &leaves {
    print_node(&mut writer, node, names)?;
  }
  print_fake_edges(&mut writer, &leaves.iter().map(|node| node.key()).collect_vec())?;

  writeln!(writer, "  }}")?;
  Ok(())
}

fn print_node<W>(mut writer: W, node: &Node, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Result<(), Report>
where
  W: Write,
{
  let key = node.key();
  let label = names[&key].clone();

  if let Some(label) = label {
    let label = label.replace('\\', "\\\\").replace('"', "\\\"");
    writeln!(writer, "    {key} [label=\"({key}) {label}\"]")?;
  } else {
    writeln!(writer, "    {key} [label=\"({key})\"]")?;
  }
  Ok(())
}

fn print_edges<W>(graph: &Graph, weights: &BTreeMap<GraphEdgeKey, Option<f64>>, mut writer: W) -> Result<(), Report>
where
  W: Write,
{
  for node in graph.get_nodes() {
    for edge_key in node.outbound() {
      let edge = graph
        .get_edge(*edge_key)
        .ok_or_else(|| make_internal_report!("Outbound edge {edge_key} not found in graph"))?;
      let source = edge.source();
      let target = edge.target();

      let weight = weights[edge_key];
      let label = weight.map(|weight| format_weight(weight, &NwkWriteOptions::default()));

      let mut attrs = Vec::new();
      if let Some(label) = label {
        attrs.push(format!("xlabel=\"{label}\""));
      }
      if let Some(weight) = weight {
        attrs.push(format!("weight=\"{weight}\""));
      }
      let attrs = attrs.join(", ");

      if !attrs.is_empty() {
        writeln!(writer, "  {source} -> {target} [{attrs}]")?;
      } else {
        writeln!(writer, "  {source} -> {target}")?;
      }
    }
  }
  Ok(())
}

fn print_fake_edges<W>(mut writer: W, node_keys: &[GraphNodeKey]) -> Result<(), Report>
where
  W: Write,
{
  let fake_edges = iproduct!(node_keys, node_keys)
    .enumerate()
    .map(|(i, (left, right))| {
      let weight = FAKE_EDGE_BASE_WEIGHT + i * FAKE_EDGE_WEIGHT_STEP;
      format!("      {left}-> {right} [style=invis, weight={weight}]")
    })
    .join("\n");

  if !fake_edges.is_empty() {
    writeln!(
      writer,
      "\n    // fake edges for alignment of nodes\n    {{\n      rank=same\n{fake_edges}\n    }}"
    )?;
  }
  Ok(())
}
