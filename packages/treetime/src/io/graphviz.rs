use crate::graph::edge::GraphEdge;
use crate::graph::graph::{Graph, SafeNode};
use crate::graph::node::{GraphNode, Node};
use crate::io::file::create_file_or_stdout;
use eyre::Report;
use itertools::{Itertools, iproduct};
use parking_lot::RwLock;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;
use std::sync::Arc;

pub fn graphviz_write_file<N, E, D>(filepath: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToGraphviz,
  E: GraphEdge + EdgeToGraphViz,
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
  E: GraphEdge + EdgeToGraphViz,
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
  E: GraphEdge + EdgeToGraphViz,
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
  print_nodes(graph, &mut writer);
  writeln!(writer)?;
  print_edges(graph, &mut writer);
  writeln!(writer, "}}")?;
  Ok(())
}

fn print_node<W, N>(mut writer: W, node: &SafeNode<N>)
where
  W: Write,
  N: GraphNode + NodeToGraphviz,
{
  let key = node.read_arc().key();
  let node = node.read_arc().payload().read_arc();
  let label = node.to_graphviz_label();

  if let Some(label) = label {
    let label = label.as_ref();
    writeln!(writer, "    {key} [label=\"({key}) {label}\"]").unwrap();
  } else {
    writeln!(writer, "    {key} [label=\"({key})\"]").unwrap();
  }
}

fn print_nodes<W, N, E, D>(graph: &Graph<N, E, D>, mut writer: W)
where
  W: Write,
  N: GraphNode + NodeToGraphviz,
  E: GraphEdge + EdgeToGraphViz,
  D: Send + Sync,
{
  writeln!(writer, "\n  subgraph roots {{").unwrap();
  let roots = graph.get_roots();
  for node in &roots {
    print_node(&mut writer, node);
  }
  print_fake_edges(&mut writer, &roots);

  writeln!(writer, "  }}\n\n  subgraph internals {{").unwrap();
  let internal = graph.get_internal_nodes();
  for node in internal {
    print_node(&mut writer, &node);
  }

  writeln!(writer, "  }}\n\n  subgraph leaves {{").unwrap();
  let leaves = graph.get_leaves();
  for node in &leaves {
    print_node(&mut writer, node);
  }
  print_fake_edges(&mut writer, &leaves);

  writeln!(writer, "  }}").unwrap();
}

fn print_edges<W, N, E, D>(graph: &Graph<N, E, D>, mut writer: W)
where
  W: Write,
  N: GraphNode + NodeToGraphviz,
  E: GraphEdge + EdgeToGraphViz,
  D: Send + Sync,
{
  graph.get_nodes().iter().for_each(&mut |node: &Arc<RwLock<Node<N>>>| {
    for edge_key in node.read().outbound() {
      let edge = graph.get_edge(*edge_key).unwrap();
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
        writeln!(writer, "  {source} -> {target} [{attrs}]").unwrap();
      } else {
        writeln!(writer, "  {source} -> {target}").unwrap();
      }
    }
  });
}

fn print_fake_edges<W, N>(mut writer: W, nodes: &[SafeNode<N>])
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
    )
    .unwrap();
  }
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

/// Defines how to display edge information when writing to GraphViz (.dot) file
pub trait EdgeToGraphViz {
  // Defines how to display label (name) of the edge in GraphViz (.dot) file
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>>;

  // Defines how to assign weight of the edge in GraphViz (.dot) file
  fn to_graphviz_weight(&self) -> Option<f64>;

  // Defines how to display additional attributes of the edge in GraphViz (.dot) file
  fn to_graphviz_attributes(&self) -> BTreeMap<String, String> {
    BTreeMap::<String, String>::new()
  }
}
