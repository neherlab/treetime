use crate::cli::treetime_cli::{MethodAnc, TreetimeAncestralArgs};
use crate::graph::graph::{Graph as GenericGraph, Weighted};
use crate::gtr::get_gtr::get_gtr;
use crate::io::fasta::read_many_fasta;
use crate::io::nwk::read_nwk;
use color_eyre::Section;
use eyre::{eyre, Report, WrapErr};
use indexmap::IndexMap;
use std::fmt::{Display, Formatter};
use std::path::Path;

pub fn run_ancestral(ancestral_args: &TreetimeAncestralArgs) -> Result<(), Report> {
  let TreetimeAncestralArgs {
    input_fastas,
    aln,
    vcf_reference,
    tree,
    gtr,
    gtr_params,
    aa,
    marginal,
    keep_overhangs,
    zero_based,
    reconstruct_tip_states,
    report_ambiguous,
    method_anc,
    outdir,
  } = ancestral_args;

  let fasta_records = read_many_fasta(input_fastas)?;

  let gtr = get_gtr(gtr)?;

  let graph = match tree {
    None => infer_graph()?,
    Some(tree) => create_graph(tree)?,
  };

  ml_anc(ancestral_args);

  Ok(())
}

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

#[derive(Clone, Debug, PartialEq)]
pub struct EdgePayload {
  weight: f64,
}

impl Weighted for EdgePayload {}

impl Display for EdgePayload {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{:1.4}", &self.weight)
  }
}

type Graph = GenericGraph<NodePayload, EdgePayload>;

fn infer_graph() -> Result<Graph, Report> {
  unimplemented!("Not implemented: Graph inference is not yet implemented. Please provide the `--tree` argument.")
}

fn create_graph(tree_path: impl AsRef<Path>) -> Result<Graph, Report> {
  let nwk_tree = read_nwk(tree_path)
    .wrap_err("When parsing input tree")
    .with_section(|| "Note: only Newick format is currently supported")?;

  let mut graph = Graph::new();

  // Insert nodes
  let mut index_map = IndexMap::<usize, usize>::new(); // Map of internal `nwk` node indices to `Graph` node indices
  for (nwk_idx, nwk_node) in nwk_tree.g.raw_nodes().iter().enumerate() {
    // Attempt to parse weight as float. If not a float, then it's a named leaf node, otherwise - internal node.
    let inserted_node_idx = match nwk_node.weight.parse::<f64>() {
      Ok(weight) => graph.add_node(NodePayload::Internal(weight)),
      Err(_) => graph.add_node(NodePayload::Leaf(nwk_node.weight.clone())),
    };

    index_map.insert(nwk_idx, inserted_node_idx);
  }

  // Insert edges
  for (nwk_idx, nwk_edge) in nwk_tree.g.raw_edges().iter().enumerate() {
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

  Ok(graph)
}

fn ml_anc(ancestral_args: &TreetimeAncestralArgs) {
  match ancestral_args.method_anc {
    MethodAnc::Ml | MethodAnc::Probabilistic => {
      if ancestral_args.marginal {
        ml_anc_marginal()
      } else {
        ml_anc_joint()
      }
    }
    MethodAnc::Parsimony | MethodAnc::Fitch => ml_anc_fitch(),
  }
}

fn ml_anc_marginal() {}

fn ml_anc_joint() {}

fn ml_anc_fitch() {}
