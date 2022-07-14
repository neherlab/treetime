use ctor::ctor;
use eyre::Report;
use itertools::Itertools;
use std::borrow::Borrow;
use std::error::Error;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::io::Write;
use std::time::Duration;
use treetime::graph::graph::Graph;
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
#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct NodePayload {
  name: String,
}

impl Display for NodePayload {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.name)
  }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct EdgePayload {
  name: String,
}

impl Display for EdgePayload {
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
  println!("{:^6} | {:^16} | {:^5}", "Node", "Parents", "Is leaf");
  graph.par_iter_breadth_first_forward(|node| {
    let is_leaf = node.is_leaf();

    let mut node_payload = node.payload_mut();

    node_payload.name = format!("{}*", &node_payload.name);
    let node_name = &node_payload.name;

    let parent_edges = node.inbound();
    let parents = parent_edges.iter().map(|parent_edge| {
      let parent_edge = parent_edge.upgrade().unwrap();
      let parent_edge_payload = parent_edge.load();
      let parent_node_payload = parent_edge.source().read().payload().clone();
      (parent_node_payload, parent_edge_payload)
    });

    let parent_names = parents.map(|(node, _)| node.name).join(", ");

    println!("{:<6} | {:<16} | {:<5}", node_name, parent_names, is_leaf);
  });

  graph.print_graph(create_file("tmp/graph2.dot")?)?;

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

  graph.build()?;

  graph.print_graph(create_file("tmp/graph.dot")?)?;

  Ok(graph)
}
