#![allow(clippy::many_single_char_names)]

use ctor::ctor;
use derive_more::Display;
use eyre::Report;
use itertools::Itertools;
use treetime::graph::graph::Graph;
use treetime::io::file::create_file;
use treetime::utils::global_init::global_init;

#[cfg(all(target_family = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
}

/// An example of node payload type
#[derive(Clone, Debug, Display)]
pub struct NodePayload {
  name: String,
}

/// An example of edge payload type
#[derive(Clone, Debug, Display)]
pub struct EdgePayload {
  name: String,
}

fn main() -> Result<(), Report> {
  let mut graph = create_example_graph()?;

  graph.write_graph(create_file("tmp/graph.dot")?)?;

  // Finalize the graph. Performs some internal bookkeeping.
  graph.build();

  println!("Traverse forward:");
  println!("{:^6} | {:^20} | {:^}", "Node", "Parents", "Is root");
  graph.iter_breadth_first_forward(|node| {
    let node_data = &node.payload.name;
    let parent_names = &node.parents.iter().map(|parent| &parent.name).join(", ");
    let is_root = &node.is_root();
    println!("{node_data:>6} | {parent_names:20} | {is_root}");
  });

  println!();

  println!("Traverse backwards:");
  println!("{:^6} | {:^20} | {:^}", "Node", "Children", "Is leaf");
  graph.iter_breadth_first_reverse(|node| {
    let node_name = &node.payload.name;
    let child_names = &node.children.iter().map(|child| &child.name).join(", ");
    let is_leaf = &node.is_leaf();
    println!("{node_name:>6} | {child_names:20} | {is_leaf}");
  });

  Ok(())
}

#[rustfmt::skip]
fn create_example_graph() -> Result<Graph<NodePayload, EdgePayload>, Report>  {
  let mut graph = Graph::<NodePayload, EdgePayload>::new();

  // Create nodes of the graph by providing their payloads. Node payload can be any arbitrary type.
  // In this case it is a `struct NodePayload` defined just above.
  //
  // At this point, nodes are not connected yet. Each insertion operation returns an index of
  // a newly created node, which can later be used for creating graph edges.
  let r0 = graph.add_node(NodePayload { name: "r0".to_owned() });
  let r1 = graph.add_node(NodePayload { name: "r1".to_owned() });
  let r2 = graph.add_node(NodePayload { name: "r2".to_owned() });
  let a = graph.add_node(NodePayload { name: "a".to_owned() });
  let b = graph.add_node(NodePayload { name: "b".to_owned() });
  let c = graph.add_node(NodePayload { name: "c".to_owned() });
  let d = graph.add_node(NodePayload { name: "d".to_owned() });
  let e = graph.add_node(NodePayload { name: "e".to_owned() });
  let f = graph.add_node(NodePayload { name: "f".to_owned() });
  let g = graph.add_node(NodePayload { name: "g".to_owned() });
  let h = graph.add_node(NodePayload { name: "h".to_owned() });
  let i = graph.add_node(NodePayload { name: "i".to_owned() });
  let j = graph.add_node(NodePayload { name: "j".to_owned() });
  let k = graph.add_node(NodePayload { name: "k".to_owned() });
  let l = graph.add_node(NodePayload { name: "l".to_owned() });
  let m = graph.add_node(NodePayload { name: "m".to_owned() });
  let n = graph.add_node(NodePayload { name: "n".to_owned() });
  let o = graph.add_node(NodePayload { name: "o".to_owned() });

  // Connect nodes pairwise. Each connection operation creates a graph edge between a pair of nodes.
  // The edge is directed from the first node to the second node, i.e. the first node is considered
  // a parent and the second node is considered a child.
  graph.add_edge(r0, r1, EdgePayload{ name: "r0->r1".to_owned() });
  graph.add_edge(r0, r2, EdgePayload{ name: "r0->r2".to_owned() });
  graph.add_edge(r1, a, EdgePayload{ name: "r1->a".to_owned() });
  graph.add_edge(r2, b, EdgePayload{ name: "r2->b".to_owned() });
  graph.add_edge(a, c, EdgePayload{ name: "a->c".to_owned() });
  graph.add_edge(a, f, EdgePayload{ name: "a->f".to_owned() });
  graph.add_edge(a, d, EdgePayload{ name: "a->d".to_owned() });
  graph.add_edge(b, d, EdgePayload{ name: "b->d".to_owned() });
  graph.add_edge(d, g, EdgePayload{ name: "d->g".to_owned() });
  graph.add_edge(b, e, EdgePayload{ name: "b->e".to_owned() });
  graph.add_edge(e, h, EdgePayload{ name: "e->h".to_owned() });
  graph.add_edge(c, g, EdgePayload{ name: "c->g".to_owned() });
  graph.add_edge(f, i, EdgePayload{ name: "f->i".to_owned() });
  graph.add_edge(f, j, EdgePayload{ name: "f->j".to_owned() });
  graph.add_edge(d, k, EdgePayload{ name: "d->k".to_owned() });
  graph.add_edge(e, l, EdgePayload{ name: "e->l".to_owned() });
  graph.add_edge(b, m, EdgePayload{ name: "b->m".to_owned() });
  graph.add_edge(r2, n, EdgePayload{ name: "r2->n".to_owned() });
  graph.add_edge(r2, o, EdgePayload{ name: "r2->o".to_owned() });
  graph.add_edge(e, o, EdgePayload{ name: "e->o".to_owned() });
  graph.add_edge(m, k, EdgePayload{ name: "m->k".to_owned() });

  graph.write_graph(create_file("tmp/graph.dot")?)?;

  Ok(graph)
}
