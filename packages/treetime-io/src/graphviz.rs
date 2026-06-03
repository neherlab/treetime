use crate::nwk::{NwkWriteOptions, format_weight};
use eyre::Report;
use itertools::{Itertools, iproduct};
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::{Graph, SafeNode};
use treetime_graph::node::GraphNode;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::make_internal_report;

pub fn graphviz_write_file<N, E, D>(filepath: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToGraphviz,
  E: GraphEdge + EdgeToGraphviz,
  D: Send + Sync,
{
  let mut f = create_file_or_stdout(filepath)?;
  graphviz_write(&mut f, graph)?;
  writeln!(f)?;
  Ok(())
}

pub fn graphviz_write_str<N, E, D>(graph: &Graph<N, E, D>) -> Result<String, Report>
where
  N: GraphNode + NodeToGraphviz,
  E: GraphEdge + EdgeToGraphviz,
  D: Send + Sync,
{
  let mut buf = Vec::new();
  graphviz_write(&mut buf, graph)?;
  Ok(String::from_utf8(buf)?)
}

pub fn graphviz_write<W, N, E, D>(mut writer: W, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  W: Write,
  N: GraphNode + NodeToGraphviz,
  E: GraphEdge + EdgeToGraphviz,
  D: Send + Sync,
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
  print_nodes(graph, &mut writer)?;
  writeln!(writer)?;
  print_edges(graph, &mut writer)?;
  writeln!(writer, "}}")?;
  Ok(())
}

fn print_node<W, N>(mut writer: W, node: &SafeNode<N>) -> Result<(), Report>
where
  W: Write,
  N: GraphNode + NodeToGraphviz,
{
  let key = node.read_arc().key();
  let node = node.read_arc().payload().read_arc();
  let label = node.to_graphviz_label();

  if let Some(label) = label {
    let label = label.as_ref();
    writeln!(writer, "    {key} [label=\"({key}) {label}\"]")?;
  } else {
    writeln!(writer, "    {key} [label=\"({key})\"]")?;
  }
  Ok(())
}

fn print_nodes<W, N, E, D>(graph: &Graph<N, E, D>, mut writer: W) -> Result<(), Report>
where
  W: Write,
  N: GraphNode + NodeToGraphviz,
  E: GraphEdge + EdgeToGraphviz,
  D: Send + Sync,
{
  writeln!(writer, "\n  subgraph roots {{")?;
  let roots = graph.get_roots();
  for node in &roots {
    print_node(&mut writer, node)?;
  }
  print_fake_edges(&mut writer, &roots)?;

  writeln!(writer, "  }}\n\n  subgraph internals {{")?;
  let internal = graph.get_internal_nodes();
  for node in internal {
    print_node(&mut writer, &node)?;
  }

  writeln!(writer, "  }}\n\n  subgraph leaves {{")?;
  let leaves = graph.get_leaves();
  for node in &leaves {
    print_node(&mut writer, node)?;
  }
  print_fake_edges(&mut writer, &leaves)?;

  writeln!(writer, "  }}")?;
  Ok(())
}

fn print_edges<W, N, E, D>(graph: &Graph<N, E, D>, mut writer: W) -> Result<(), Report>
where
  W: Write,
  N: GraphNode + NodeToGraphviz,
  E: GraphEdge + EdgeToGraphviz,
  D: Send + Sync,
{
  for node in graph.get_nodes() {
    for edge_key in node.read().outbound() {
      let edge = graph
        .get_edge(*edge_key)
        .ok_or_else(|| make_internal_report!("Outbound edge {edge_key} not found in graph"))?;
      let source = edge.read_arc().source();
      let target = edge.read_arc().target();
      let payload = edge.read_arc().payload().read_arc();

      let label = payload.to_graphviz_label();
      let weight = payload.to_graphviz_weight();

      let mut attrs = Vec::new();
      if let Some(label) = label {
        let label = label.as_ref();
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

fn print_fake_edges<W, N>(mut writer: W, nodes: &[SafeNode<N>]) -> Result<(), Report>
where
  W: Write,
  N: GraphNode + NodeToGraphviz,
{
  // Fake edges needed to align a set of nodes beautifully
  let node_keys = nodes.iter().map(|node| node.read().key()).collect_vec();

  let fake_edges = iproduct!(&node_keys, &node_keys)
    .enumerate()
    .map(|(i, (left, right))| {
      let weight = 1000 + i * 100;
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

/// Defines how to display node information when writing to GraphViz (.dot) file
pub trait NodeToGraphviz {
  // Defines how to display label (name) of the node in GraphViz (.dot) file
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>>;

  // Defines how to display additional attributes of the node in GraphViz (.dot) file
  fn to_graphviz_attributes(&self) -> BTreeMap<String, String> {
    BTreeMap::<String, String>::new()
  }
}

/// Defines how to display edge information when writing to GraphViz (.dot) file.
/// Default implementations derive label and weight from `HasBranchLength`.
pub trait EdgeToGraphviz: HasBranchLength {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .branch_length()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.branch_length()
  }

  fn to_graphviz_attributes(&self) -> BTreeMap<String, String> {
    BTreeMap::<String, String>::new()
  }
}
