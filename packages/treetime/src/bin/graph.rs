use ctor::ctor;
use eyre::Report;
use itertools::Itertools;
use treetime::graph::graph::Graph;
use treetime::utils::global_init::global_init;

#[cfg(all(target_family = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
}

/// An example of node payload type
#[derive(Clone, Debug)]
pub struct NodePayload {
  name: &'static str,
}

fn main() -> Result<(), Report> {
  let mut graph = Graph::new();

  // Create nodes of the graph by providing their payloads. Node payload can be any arbitrary type.
  // In this case it is a `struct NodePayload` defined just above.
  //
  // At this point, nodes are not connected yet. Each insertion operation returns an index of
  // a newly created node, which can later be used for creating graph edges.
  let root1 = graph.insert(NodePayload { name: "root1" });
  let root2 = graph.insert(NodePayload { name: "root2" });
  let a = graph.insert(NodePayload { name: "a" });
  let b = graph.insert(NodePayload { name: "b" });
  let c = graph.insert(NodePayload { name: "c" });
  let d = graph.insert(NodePayload { name: "d" });
  let e = graph.insert(NodePayload { name: "e" });
  let f = graph.insert(NodePayload { name: "f" });
  let g = graph.insert(NodePayload { name: "g" });
  let h = graph.insert(NodePayload { name: "h" });

  // Connect nodes pairwise. Each connection operation creates a graph edge between a pair of nodes.
  // The edge is directed from the first node to the second node, i.e. the first node is considered
  // a parent and the second node is considered a child.
  graph.connect(root1, a);
  graph.connect(root2, b);
  graph.connect(a, c);
  graph.connect(a, f);
  graph.connect(a, d);
  graph.connect(b, d);
  graph.connect(d, g);
  graph.connect(b, e);
  graph.connect(e, h);

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
  println!("{:^6} | {:^20} | {:^}", "Node", "Parents", "Is leaf");
  graph.iter_breadth_first_reverse(|node| {
    let node_name = &node.payload.name;
    let child_names = &node.children.iter().map(|child| &child.name).join(", ");
    let is_leaf = &node.is_leaf();
    println!("{node_name:>6} | {child_names:20} | {is_leaf}");
  });

  Ok(())
}
