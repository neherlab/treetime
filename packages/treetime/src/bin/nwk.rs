use bio::io::newick;
use clap::{AppSettings, Parser};
use ctor::ctor;
use eyre::{eyre, Report};
use indexmap::IndexMap;
use std::borrow::Borrow;
use std::fmt::{Debug, Display, Formatter};
use std::hash::Hash;
use std::path::PathBuf;
use treetime::graph::graph::Graph;
use treetime::io::file::create_file;
use treetime::io::fs::read_file_to_string;
use treetime::utils::global_init::global_init;

#[cfg(all(target_os = "linux", target_arch = "x86_64"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[ctor]
fn init() {
  global_init();
}

/// An example of node payload type
#[derive(Clone, Debug, PartialEq)]
pub enum NodePayload {
  Leaf(String),
  Internal(f64),
}

impl Display for NodePayload {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    match self {
      NodePayload::Leaf(id) => write!(f, "{id}"),
      NodePayload::Internal(weight) => write!(f, "{weight:1.4}"),
    }
  }
}

/// An example of edge payload type
#[derive(Clone, Debug, PartialEq)]
pub struct EdgePayload {
  weight: f64,
}

impl Display for EdgePayload {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{:1.4}", &self.weight)
  }
}

#[derive(Parser, Debug)]
#[clap(name = "treetime-nwk", trailing_var_arg = true)]
#[clap(author, version)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(verbatim_doc_comment)]
pub struct TreetimeNwkArgs {
  #[clap(long, short = 'i')]
  input_nwk: PathBuf,

  #[clap(long, short = 'o')]
  output_dot: PathBuf,
}

fn main() -> Result<(), Report> {
  let TreetimeNwkArgs { input_nwk, output_dot } = TreetimeNwkArgs::parse();

  let mut graph = Graph::<NodePayload, EdgePayload>::new();

  let nwk_str = read_file_to_string(&input_nwk)?;
  let tree = newick::read(nwk_str.as_bytes())?;

  // Insert nodes
  let mut index_map = IndexMap::<usize, usize>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  for (nwk_idx, nwk_node) in tree.g.raw_nodes().iter().enumerate() {
    // Attempt to parse weight as float. If not a float, then it's a named leaf node, otherwise - internal node.
    let inserted_node_idx = match nwk_node.weight.parse::<f64>() {
      Ok(weight) => graph.add_node(NodePayload::Internal(weight)),
      Err(_) => graph.add_node(NodePayload::Leaf(nwk_node.weight.clone())),
    };

    index_map.insert(nwk_idx, inserted_node_idx);
  }

  // Insert edges
  for (nwk_idx, nwk_edge) in tree.g.raw_edges().iter().enumerate() {
    let weight: f64 = nwk_edge.weight as f64;
    let source: usize = nwk_edge.source().index();
    let target: usize = nwk_edge.target().index();

    let source = index_map
      .get(&source)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {source} not found."))?;

    let target = index_map
      .get(&target)
      .ok_or_else(|| eyre!("When inserting edge {nwk_idx}: Node with index {target} not found."))?;

    graph.add_edge(*source, *target, EdgePayload { weight });
  }

  graph.print_graph(create_file(&output_dot)?)?;

  Ok(())
}
