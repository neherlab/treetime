#![allow(clippy::many_single_char_names)]

use ctor::ctor;
use eyre::Report;
use itertools::Itertools;
use std::borrow::Borrow;
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::io::Write;
use std::time::Duration;
use treetime::graph::breadth_first::GraphTraversalContinuation;
use treetime::graph::edge::{GraphEdge, Weighted};
use treetime::graph::graph::{Graph, NodeEdgePair};
use treetime::graph::node::{GraphNode, Named};
use treetime::io::file::create_file;
use treetime::utils::global_init::global_init;

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
}

/// An example of node payload type
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Node {
  name: String,
}

impl Named for Node {
  fn name(&self) -> &str {
    &self.name
  }
}

impl GraphNode for Node {}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.name)
  }
}

/// An example of edge payload type
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Edge {
  name: String,
}

impl Weighted for Edge {
  fn weight(&self) -> f64 {
    1.0
  }
}

impl GraphEdge for Edge {}

impl Display for Edge {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", &self.name)
  }
}

fn main() -> Result<(), Report> {
  let mut graph = create_example_graph()?;

  #[cfg(debug_assertions)]
  {
    use parking_lot::deadlock;
    use std::thread;
    use std::time::Duration;

    // Create a background thread which checks for deadlocks every 10s
    thread::spawn(move || loop {
      thread::sleep(Duration::from_secs(5));
      let deadlocks = deadlock::check_deadlock();
      println!("{} deadlocks detected", deadlocks.len());

      if deadlocks.is_empty() {
        continue;
      }

      for (i, threads) in deadlocks.iter().enumerate() {
        println!("Deadlock #{}", i);
        for t in threads {
          println!("Thread Id {:#?}", t.thread_id());
          println!("{:#?}", t.backtrace());
        }
      }
    });
  }

  println!("Traverse forward:");
  println!(
    "{:^6} | {:^16} | {:^7} | {:^7}",
    "Node", "Parents", "Is leaf", "Is root"
  );
  graph.par_iter_breadth_first_forward(|node| {
    // Mutable access to node payload
    node.payload.name = format!("{}*", &node.payload.name);

    // Read-write access to list pairs of `(edge, parent)`, where `edge is the edge leading to the `parent`.
    // Note: order of writes is not specified. Write access is exclusive and blocks other threads.
    // Note: there is no access to children, because they are not guaranteed to be processed at this point yet.
    let parent_names = node
      .parents
      .iter()
      .map(|NodeEdgePair { edge, node: parent }| {
        let parent = parent.read();
        parent.name.clone()
      })
      .join(", ");

    println!(
      "{:<6} | {:<16} | {:<7} | {:<7}",
      node.payload.name, parent_names, node.is_leaf, node.is_root
    );

    GraphTraversalContinuation::Continue
  });

  graph.print_graph(create_file("tmp/graph2.dot")?)?;

  println!();

  println!("Traverse backward:");
  println!(
    "{:^6} | {:^16} | {:^7} | {:^7}",
    "Node", "Children", "Is leaf", "Is root"
  );
  graph.par_iter_breadth_first_backward(|node| {
    // Mutable access to node payload
    node.payload.name = format!("{}*", &node.payload.name);

    // Access to list pairs of `(edge, child)`, where `edge` is the edge leading to that `child`
    // Note: order of writes is not specified. Write access is exclusive and blocks other threads.
    // Note: there is no access to parents, because they are not guaranteed to be processed at this point yet.
    let child_names = node
      .children
      .iter()
      .map(|NodeEdgePair { edge, node: child }| {
        let child = child.read();
        child.name.clone()
      })
      .join(", ");

    println!(
      "{:<6} | {:<16} | {:<7} | {:<7}",
      node.payload.name, child_names, node.is_leaf, node.is_root
    );

    GraphTraversalContinuation::Continue
  });

  Ok(())
}

#[rustfmt::skip]
fn create_example_graph() -> Result<Graph<Node, Edge>, Report>  {
  let mut graph = Graph::<Node, Edge>::new();

  // Create nodes of the graph by providing their payloads. Node payload can be any arbitrary type.
  // In this case it is a `struct NodePayload` defined just above.
  //
  // At this point, nodes are not connected yet. Each insertion operation returns an index of
  // a newly created node, which can later be used for creating graph edges.
  let r1 = graph.add_node(Node { name: "r1".to_owned() });
  let r2 = graph.add_node(Node { name: "r2".to_owned() });
  let  a = graph.add_node(Node { name:  "a".to_owned() });
  let  b = graph.add_node(Node { name:  "b".to_owned() });
  let  c = graph.add_node(Node { name:  "c".to_owned() });
  let  d = graph.add_node(Node { name:  "d".to_owned() });
  let  e = graph.add_node(Node { name:  "e".to_owned() });
  let  f = graph.add_node(Node { name:  "f".to_owned() });
  let  g = graph.add_node(Node { name:  "g".to_owned() });
  let  h = graph.add_node(Node { name:  "h".to_owned() });
  let  i = graph.add_node(Node { name:  "i".to_owned() });
  let  j = graph.add_node(Node { name:  "j".to_owned() });
  let  k = graph.add_node(Node { name:  "k".to_owned() });
  let  l = graph.add_node(Node { name:  "l".to_owned() });
  let  m = graph.add_node(Node { name:  "m".to_owned() });
  let  n = graph.add_node(Node { name:  "n".to_owned() });
  let  o = graph.add_node(Node { name:  "o".to_owned() });

  // Connect nodes pairwise. Each connection operation creates a graph edge between a pair of nodes.
  // The edge is directed from the first node to the second node, i.e. the first node is considered
  // a parent and the second node is considered a child.
  graph.add_edge(r1, a, Edge { name: "r1->a".to_owned() });
  graph.add_edge(r2, b, Edge { name: "r2->b".to_owned() });
  graph.add_edge(a,  c, Edge { name:  "a->c".to_owned() });
  graph.add_edge(a,  f, Edge { name:  "a->f".to_owned() });
  graph.add_edge(a,  d, Edge { name:  "a->d".to_owned() });
  graph.add_edge(b,  d, Edge { name:  "b->d".to_owned() });
  graph.add_edge(d,  g, Edge { name:  "d->g".to_owned() });
  graph.add_edge(b,  e, Edge { name:  "b->e".to_owned() });
  graph.add_edge(e,  h, Edge { name:  "e->h".to_owned() });
  graph.add_edge(c,  g, Edge { name:  "c->g".to_owned() });
  graph.add_edge(f,  i, Edge { name:  "f->i".to_owned() });
  graph.add_edge(f,  j, Edge { name:  "f->j".to_owned() });
  graph.add_edge(d,  k, Edge { name:  "d->k".to_owned() });
  graph.add_edge(e,  l, Edge { name:  "e->l".to_owned() });
  graph.add_edge(b,  m, Edge { name:  "b->m".to_owned() });
  graph.add_edge(r2, n, Edge { name: "r2->n".to_owned() });
  graph.add_edge(r2, o, Edge { name: "r2->o".to_owned() });
  graph.add_edge(e,  o, Edge { name:  "e->o".to_owned() });
  graph.add_edge(m,  k, Edge { name:  "m->k".to_owned() });

  graph.build()?;

  graph.print_graph(create_file("tmp/graph.dot")?)?;

  Ok(graph)
}
