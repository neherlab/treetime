use crate::branch_lengths::branch_lengths_or_zero;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::{GtrModel, GtrModelName, JC69Params};
use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
use crate::partition::create::build_sparse_reconstruction;
use crate::partition::inferred_rate::InferredRate;
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::resolve::PartitionLayout;
use crate::partition::result::{
  AlignmentResult, AlignmentSequences, PartitionOutcome, PartitionResult, alignment_results,
};
use crate::partition::set::layouts_or_whole;
use crate::partition::split::{AlignmentInfo, PartitionData, SplitAlignments, split_alignments};
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
  if needs_sequences && input.alignments.is_empty() {
    return Err(OperationError::InvalidInput(make_report!(
      "Sequences required for --collapse-empty or --merge-shared-mutations"
    )));
  }
  let SplitAlignments { alignments, partitions } =
    split_alignments(take(&mut input.alignments), take(&mut input.layouts)).map_err(OperationError::classify)?;
  let fitch_branch_lengths = branch_lengths_or_zero(&branch_lengths);
  let rate = InferredRate::for_partition_count(partitions.len());
  let mut partitions = partitions
    .into_iter()
    .map(
      |PartitionData {
         layout,
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
        Ok(reconstruction.with_layout(layout))
      },
    )
    .collect::<Result<Vec<_>, Report>>()
    .map_err(OperationError::classify)?;

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

  let (alignments, partitions) =
    prune_results(&input.graph, alignments, partitions).map_err(OperationError::InferenceFailed)?;

  Ok(PruneOutput {
    graph: input.graph,
    alignments,
    partitions,
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
  pub alignments: Vec<SequenceSource>,
  pub layouts: Vec<PartitionLayout>,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug)]
pub struct PruneOutput {
  pub graph: Graph,
  pub alignments: Vec<AlignmentResult>,
  pub partitions: Vec<PartitionResult>,
  pub names: NodeNames,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

fn prune_results(
  graph: &Graph,
  alignments: Vec<AlignmentInfo>,
  partitions: Vec<SparseReconstruction>,
) -> Result<(Vec<AlignmentResult>, Vec<PartitionResult>), Report> {
  let (results, outcomes): (Vec<_>, Vec<_>) = layouts_or_whole(
    partitions
      .iter()
      .map(|partition| (partition.layout(), partition.partition.length)),
  )
  .into_iter()
  .zip(partitions)
  .map(|(layout, SparseReconstruction { partition, gtr, .. })| {
    let outcome = PartitionOutcome {
      alignment: layout.alignment,
      indel_rate: None,
      sequences: Some(fitch_sequences(graph, &partition)?),
    };
    Ok((
      PartitionResult::new(&layout, GtrModelName::JC69, Some(gtr), None),
      outcome,
    ))
  })
  .collect::<Result<Vec<_>, Report>>()?
  .into_iter()
  .unzip();
  Ok((alignment_results(alignments, outcomes)?, results))
}

fn fitch_sequences(graph: &Graph, partition: &PartitionMarginalSparse) -> Result<AlignmentSequences, Report> {
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
  AlignmentSequences::new(mutations, &partition.alphabet)
}
