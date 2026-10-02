use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
use crate::partition::create::build_sparse_partition;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::progress::LogSink;
use crate::prune::prune::prune_nodes;
use crate::seq::alignment::node_seq_inputs;
use eyre::eyre;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use std::mem::take;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AlignmentRecord;

pub fn run(
  params: &PruneParams,
  mut input: PruneInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  cancel: &dyn Cancel,
  log: &dyn LogSink,
) -> Result<PruneOutput, OperationError> {
  cancel.check()?;

  let mut names = names.clone();
  let mut branch_lengths = take(&mut input.branch_lengths);

  let needs_sequences = params.prune_empty || params.merge_shared_mutations;
  let (mut partitions, gtr) = if needs_sequences {
    let sequences = take(&mut input.sequences).ok_or_else(|| {
      OperationError::InvalidInput(eyre!(
        "Sequences required for --prune-empty or --merge-shared-mutations"
      ))
    })?;
    let node_inputs = node_seq_inputs(&input.graph, &names, sequences);
    let (partition, _, gtr) = build_sparse_partition(
      GtrModelName::JC69,
      &input.graph,
      0,
      input.alphabet,
      &node_inputs,
      &branch_lengths_or_zero(&branch_lengths),
      log,
    )?;
    (vec![partition], Some(gtr))
  } else {
    (vec![], None)
  };

  prune_nodes(
    &mut input.graph,
    &mut partitions,
    params.prune_short,
    params.prune_empty,
    &params.node_names,
    &names,
    &mut branch_lengths,
  )?;

  if params.merge_shared_mutations {
    merge_shared_mutation_branches(&mut input.graph, &mut partitions, &mut branch_lengths)?;
    input.graph.build()?;
    names = assign_node_names(names, &input.graph)?;
  }

  Ok(PruneOutput {
    graph: input.graph,
    gtr,
    partitions,
    names,
    branch_lengths,
  })
}

pub struct PruneParams {
  pub prune_short: Option<f64>,
  pub prune_empty: bool,
  pub merge_shared_mutations: bool,
  pub node_names: BTreeSet<String>,
}

pub struct PruneInput {
  pub graph: Graph,
  pub alphabet: Alphabet,
  pub sequences: Option<Vec<AlignmentRecord>>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug, Serialize)]
pub struct PruneOutput {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub partitions: Vec<PartitionMarginalSparse>,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}
