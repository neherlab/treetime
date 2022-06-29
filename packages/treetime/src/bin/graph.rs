use ctor::ctor;
use eyre::Report;
use treetime::graph::graph::{Graph, GraphNode};
use treetime::utils::global_init::global_init;

#[cfg(all(target_family = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
}

#[derive(Clone, Debug)]
pub struct NodePayload {
  data: &'static str,
}

fn main() -> Result<(), Report> {
  let mut graph = Graph::new();
  let root = graph.insert(NodePayload { data: "root" });
  let a = graph.insert(NodePayload { data: "a" });
  let b = graph.insert(NodePayload { data: "b" });
  let c = graph.insert(NodePayload { data: "c" });
  let d = graph.insert(NodePayload { data: "d" });
  let e = graph.insert(NodePayload { data: "e" });
  let f = graph.insert(NodePayload { data: "f" });
  let g = graph.insert(NodePayload { data: "g" });
  let h = graph.insert(NodePayload { data: "h" });

  graph.connect(root, a);
  graph.connect(root, b);

  graph.connect(a, c);
  graph.connect(a, f);

  graph.connect(a, d);
  graph.connect(b, d);
  graph.connect(d, g);

  graph.connect(b, e);
  graph.connect(e, h);

  for node in graph.iter_breadth_first() {
    println!("{}", &node.payload.data);
  }

  Ok(())
}
