use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::{GtrModel, GtrModelName, JC69Params};
use crate::gtr::gtr::GTR;
use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
use crate::partition::create::build_sparse_reconstruction;
use crate::partition::inferred_rate::InferredRate;
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::result::{PartitionResult, PartitionSequences, PartitionSource, SplitSource, split_sources};
use crate::progress::LogSink;
use crate::prune::prune::prune_nodes;
use crate::seq::alignment::SequenceSource;
use crate::seq::mutation::{MutationTrack, SequenceMutations};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use std::mem::take;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_utils::make_report;

pub fn run(
  params: &PruneParams,
  mut input: PruneInput,
  cancel: &dyn Cancel,
  log: &dyn LogSink,
) -> Result<PruneOutput, OperationError> {
  cancel.check().map_err(OperationError::classify)?;

  let mut branch_lengths = take(&mut input.branch_lengths);

  let needs_sequences = params.prune_empty || params.merge_shared_mutations;
  if needs_sequences && input.sources.is_empty() {
    return Err(OperationError::InvalidInput(make_report!(
      "Sequences required for --prune-empty or --merge-shared-mutations"
    )));
  }
  let sources = if needs_sequences {
    split_sources(take(&mut input.sources))
  } else {
    vec![]
  };
  let fitch_branch_lengths = branch_lengths_or_zero(&branch_lengths);
  let rate = InferredRate::for_file_count(sources.len());
  let (sources, mut partitions): (Vec<_>, Vec<_>) = sources
    .into_iter()
    .map(
      |SplitSource {
         source,
         alphabet,
         nodes,
       }| {
        let reconstruction = build_sparse_reconstruction(
          &GtrModel::JC69(JC69Params::default()),
          rate,
          &input.graph,
          alphabet,
          nodes,
          &fitch_branch_lengths,
          log,
        )?;
        Ok((source, reconstruction))
      },
    )
    .collect::<Result<Vec<_>, Report>>()
    .map_err(OperationError::classify)?
    .into_iter()
    .unzip();

  prune_nodes(
    &mut input.graph,
    &mut partitions,
    params.prune_short,
    params.prune_empty,
    &params.node_keys,
    &mut branch_lengths,
  )
  .map_err(OperationError::classify)?;

  if params.merge_shared_mutations {
    merge_shared_mutation_branches(&mut input.graph, &mut input.names, &mut partitions, &mut branch_lengths)
      .map_err(OperationError::classify)?;
    input.graph.build().map_err(OperationError::classify)?;
  }

  let results = sources
    .into_iter()
    .zip(partitions)
    .map(|(source, SparseReconstruction { partition, gtr, .. })| prune_result(&input.graph, source, &partition, gtr))
    .collect::<Result<Vec<_>, Report>>()
    .map_err(OperationError::InferenceFailed)?;

  Ok(PruneOutput {
    graph: input.graph,
    results,
    names: input.names,
    branch_lengths,
  })
}

pub struct PruneParams {
  pub prune_short: Option<f64>,
  pub prune_empty: bool,
  pub merge_shared_mutations: bool,
  pub node_keys: BTreeSet<GraphNodeKey>,
}

pub struct PruneInput {
  pub graph: Graph,
  pub names: NodeNames,
  pub sources: Vec<SequenceSource>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug)]
pub struct PruneOutput {
  pub graph: Graph,
  pub results: Vec<PartitionResult>,
  pub names: NodeNames,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

fn prune_result(
  graph: &Graph,
  source: PartitionSource,
  partition: &PartitionMarginalSparse,
  gtr: GTR,
) -> Result<PartitionResult, Report> {
  let edge_mutations = graph
    .get_edges()
    .map(|edge| {
      let edge_key = edge.key();
      Ok((
        edge_key,
        partition.edge_fitch_mutations(edge_key, &MutationTrack::Nucleotide)?,
      ))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  let mutations = SequenceMutations {
    root_sequence: partition.fitch_root_sequence(),
    edge_mutations,
  };
  Ok(PartitionResult {
    number: source.number,
    path: source.path,
    alphabet: partition.alphabet.clone(),
    length: source.length,
    mask: source.mask,
    model_name: GtrModelName::JC69,
    gtr: Some(gtr),
    indel_rate: None,
    log_lh: None,
    sequences: Some(PartitionSequences::new(mutations, &partition.alphabet)?),
  })
}
