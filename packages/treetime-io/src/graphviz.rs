use crate::nwk::{NwkWriteOptions, format_weight};
use eyre::{Report, WrapErr};
use itertools::{Itertools, iproduct};
use std::collections::BTreeMap;
use std::fmt::Write;
use std::io::Write as _;
use std::path::Path;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNodeKey, Node};
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::make_internal_report;

const FAKE_EDGE_BASE_WEIGHT: usize = 1000;
const FAKE_EDGE_WEIGHT_STEP: usize = 100;

pub fn graphviz_write_file(
  filepath: impl AsRef<Path>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let text = graphviz_write_str(graph, names, weights)?;
  let mut f = create_file_or_stdout(filepath)?;
  writeln!(f, "{text}").wrap_err_with(|| format!("When writing Graphviz file '{}'", filepath.display()))
}

fn graphviz_write_str(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<String, Report> {
  let mut text = String::new();
  graphviz_write(&mut text, graph, names, weights)?;
  Ok(text)
}

fn graphviz_write<W>(
  mut writer: W,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  weights: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<(), Report>
where
  W: Write,
{
  write!(
    writer,
    r#"
digraph Phylogeny {{
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
  for node in graph.get_internal_nodes() {
    print_node(&mut writer, node, names)?;
  }

  writeln!(writer, "  }}\n\n  subgraph leaves {{")?;
  let leaves = graph.get_leaves().collect::<Vec<_>>();
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
