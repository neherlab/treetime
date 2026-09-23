use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::ancestral::pipeline::SparseReconstruction;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr};
use crate::gtr::gtr::GTR;
use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::prune::prune::prune_nodes;
use crate::seq::alignment::node_seq_inputs;
use itertools::{Itertools, izip};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
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
) -> Result<PruneOutput, OperationError> {
  cancel.check()?;

  let mut names = names.clone();
  let mut branch_lengths = std::mem::take(&mut input.branch_lengths);

  let needs_sequences = params.prune_empty || params.merge_shared_mutations;
  let (mut partitions, gtrs, node_states): (Vec<_>, Vec<_>, Vec<_>) = if needs_sequences {
    let sequences = std::mem::take(&mut input.sequences).ok_or_else(|| {
      OperationError::InvalidInput(eyre::eyre!(
        "Sequences required for --prune-empty or --merge-shared-mutations"
      ))
    })?;
    let node_inputs = node_seq_inputs(&input.graph, &names, sequences);
    let created = create_marginal_partition(
      &input.graph,
      0,
      input.alphabet.clone(),
      &node_inputs,
      GtrModelName::JC69,
      None,
      &branch_lengths_or_zero(&branch_lengths),
    )?;
    let (partition, gtr, node_states) = match created.partition {
      MarginalPartition::Sparse(partition, node_states) => (partition, created.gtr, node_states),
      MarginalPartition::Dense(_) => {
        let gtr = get_gtr_by_name(GtrModelName::JC69)?;
        log_gtr(&gtr, GtrModelName::JC69);
        let fitch =
          crate::ancestral::fitch::create_fitch_partition(&input.graph, 0, input.alphabet.clone(), &node_inputs)?;
        let (partition, node_states) = fitch.into_marginal_sparse(&input.graph)?;
        (partition, gtr, node_states)
      },
    };
    (vec![partition], vec![gtr], vec![node_states])
  } else {
    (vec![], vec![], vec![])
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

  let gtr = gtrs.first().cloned();
  let partitions = izip!(partitions, gtrs, node_states)
    .map(|(partition, gtr, node_states)| SparseReconstruction::seeded(partition, gtr, node_states))
    .collect_vec();

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
  pub partitions: Vec<SparseReconstruction>,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}
