use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::ancestral::params::MethodAncestral;
use crate::ancestral::plan::{ReconOptions, ReconstructedPartition, reconstruct_partition, resolve_plan};
use crate::ancestral::reconstruction::ReconstructedSequences;
use crate::ancestral::sample::SampleMode;
use crate::ancestral::tip_states::TipStates;
use crate::cancel::Cancel;
use crate::error::OperationError;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::marginal::dense::partition::{DenseMarginalEdges, PartitionMarginalDense};
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::partition::marginal::sparse::partition::{PartitionMarginalSparse, SparseMarginalEdges};
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::DenseNodeState;
use crate::partition::storage::sparse::SparseNodeState;
use crate::progress::{LogSink, StageSink};
use crate::seq::alignment::AncestralInput;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use eyre::Report;
use rand::RngCore;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::AsciiChar;
use treetime_primitives::LogLh;
use treetime_primitives::Seq;
use treetime_utils::sync::random::get_random_number_generator;

pub fn run(
  params: &AncestralParams,
  input: &AncestralInput,
  alphabet: Alphabet,
  mask: Vec<bool>,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<AncestralOutputFull, OperationError> {
  let plan = resolve_plan(params)?;
  let options = ReconOptions::new(
    params.include_leaves,
    params.impute_missing_data,
    params.sample_from_profile,
  );
  let mut rng = get_random_number_generator(params.seed);
  let branch_lengths = branch_lengths_or_zero(&input.branch_lengths());
  let ReconstructedPartition {
    partition,
    emitted_nodes,
  } = reconstruct_partition(
    &input.graph,
    &plan,
    0,
    alphabet,
    &input.nodes,
    &branch_lengths,
    &options,
    &mut rng,
    cancel,
    stages,
    log,
  )?;
  Ok(AncestralOutputFull {
    output: AncestralOutput {
      gtr: partition.gtr().cloned(),
      model_name: params.model,
      mask,
      emitted_nodes,
    },
    partition: Some(partition),
  })
}

pub struct AncestralParams {
  pub method: MethodAncestral,
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub include_leaves: bool,
  pub impute_missing_data: bool,
  pub gtr_iterations: usize,
  pub site_specific_gtr: bool,
  pub seed: Option<u64>,
  pub sample_from_profile: SampleMode,
}

pub struct AncestralOutputFull {
  pub output: AncestralOutput,
  pub partition: Option<AncestralPartition>,
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum AncestralPartition {
  Fitch(PartitionFitch),
  Sparse {
    family: SparseReconstruction,
    sequences: BTreeMap<GraphNodeKey, Seq>,
  },
  Dense {
    family: DenseReconstruction,
    sequences: BTreeMap<GraphNodeKey, Seq>,
  },
}

impl AncestralPartition {
  pub(crate) fn gtr(&self) -> Option<&GTR> {
    match self {
      Self::Fitch(_) => None,
      Self::Sparse { family, .. } => Some(&family.gtr),
      Self::Dense { family, .. } => Some(&family.gtr),
    }
  }

  pub fn sequence_length(&self) -> usize {
    match self {
      Self::Fitch(partition) => partition.sequence_length(),
      Self::Sparse { family, .. } => family.sequence_length(),
      Self::Dense { family, .. } => family.sequence_length(),
    }
  }

  pub fn ambiguous_char(&self) -> AsciiChar {
    match self {
      Self::Fitch(partition) => partition.ambiguous_char(),
      Self::Sparse { family, .. } => family.ambiguous_char(),
      Self::Dense { family, .. } => family.ambiguous_char(),
    }
  }

  pub fn augur_node_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Fitch(partition) => partition.node_sequence(node_key),
      Self::Sparse { sequences, .. } | Self::Dense { sequences, .. } => sequences[&node_key].clone(),
    }
  }

  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    match self {
      Self::Fitch(partition) => partition.root_sequence(graph),
      Self::Sparse { family, .. } => family.root_sequence(graph),
      Self::Dense { family, .. } => family.root_sequence(graph),
    }
  }

  pub fn augur_root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    Ok(self.augur_node_sequence(graph.root_key()?))
  }

  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Fitch(partition) => partition.edge_subs(graph, edge_key),
      Self::Sparse { family, .. } => family.edge_subs(edge_key),
      Self::Dense { family, .. } => family.edge_subs(graph, edge_key),
    }
  }

  pub(crate) fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Fitch(partition) => partition.edge_indels(edge_key),
      Self::Sparse { family, .. } => family.edge_indels(edge_key),
      Self::Dense { family, .. } => family.edge_indels(edge_key),
    }
  }

  pub fn edge_mutations(
    &self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: &MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(graph, edge_key)?, &self.edge_indels(edge_key), track)
  }
}

#[derive(Clone, Debug, Serialize)]
pub struct SparseReconstruction {
  pub partition: PartitionMarginalSparse,
  pub gtr: GTR,
  pub node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  pub edges: SparseMarginalEdges,
}

impl SparseReconstruction {
  pub fn seeded(
    partition: PartitionMarginalSparse,
    gtr: GTR,
    node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  ) -> Self {
    Self {
      partition,
      gtr,
      node_states,
      edges: SparseMarginalEdges::default(),
    }
  }

  pub(crate) fn sequence_length(&self) -> usize {
    self.partition.length
  }

  pub fn edge_subs(&self, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    self.partition.edge_subs(&self.edges.estimates, edge_key)
  }

  pub(crate) fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    self.partition.edge_effective_length(graph, edge_key)
  }

  pub(crate) fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    self
      .partition
      .create_edge_contribution(&self.gtr, &self.edges.backward, &self.edges.forward, edge_key)
  }

  pub(crate) fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.partition.edge_indel_count(edge_key)
  }

  pub fn root_sequence(&self, _graph: &Graph) -> Result<Seq, Report> {
    Ok(self.partition.root_sequence())
  }

  pub(crate) fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    self.partition.edge_indels(edge_key)
  }

  pub fn edge_mutations(&self, edge_key: GraphEdgeKey, track: &MutationTrack) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(edge_key)?, &self.edge_indels(edge_key), track)
  }

  pub(crate) fn reconstruct_sequences(
    &self,
    graph: &Graph,
    tips: TipStates,
    sample_mode: SampleMode,
    rng: &mut dyn RngCore,
  ) -> Result<ReconstructedSequences, Report> {
    self
      .partition
      .reconstruct_sequences(graph, &self.node_states, &self.edges.forward, tips, sample_mode, rng)
  }

  fn ambiguous_char(&self) -> AsciiChar {
    self.partition.alphabet.unknown()
  }

  pub fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    let Self {
      partition,
      gtr,
      node_states,
      ..
    } = self;
    let MarginalUpdate {
      node_states,
      edges,
      log_lh,
    } = partition.marginal_update(&gtr, graph, branch_lengths, &node_states)?;
    Ok((
      Self {
        partition,
        gtr,
        node_states,
        edges,
      },
      log_lh,
    ))
  }
}

#[derive(Clone, Debug, Serialize)]
pub struct DenseReconstruction {
  pub(crate) partition: PartitionMarginalDense,
  pub(crate) gtr: GTR,
  pub(crate) node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  pub(crate) edges: DenseMarginalEdges,
}

impl DenseReconstruction {
  pub(crate) fn seeded(partition: PartitionMarginalDense, gtr: GTR) -> Self {
    Self {
      partition,
      gtr,
      node_states: BTreeMap::new(),
      edges: DenseMarginalEdges::default(),
    }
  }

  pub(crate) fn sequence_length(&self) -> usize {
    self.partition.length
  }

  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    self.partition.edge_subs(&self.node_states, graph, edge_key)
  }

  pub(crate) fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    self.partition.edge_effective_length(&self.node_states, graph, edge_key)
  }

  pub(crate) fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> OptimizationContribution {
    self
      .partition
      .create_edge_contribution(&self.gtr, &self.edges.backward, &self.edges.forward, edge_key)
  }

  pub(crate) fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    self.partition.edge_indel_count(&self.edges.estimates, edge_key)
  }

  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    self.partition.root_sequence(&self.node_states, graph)
  }

  pub(crate) fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    self.partition.edge_indels(&self.edges.estimates, edge_key)
  }

  pub fn edge_mutations(
    &self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: &MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(graph, edge_key)?, &self.edge_indels(edge_key), track)
  }

  pub(crate) fn reconstruct_sequences(
    &self,
    graph: &Graph,
    tips: TipStates,
    sample_mode: SampleMode,
    rng: &mut dyn RngCore,
  ) -> Result<ReconstructedSequences, Report> {
    self
      .partition
      .reconstruct_sequences(graph, &self.node_states, tips, sample_mode, rng)
  }

  fn ambiguous_char(&self) -> AsciiChar {
    self.partition.alphabet.unknown()
  }

  pub(crate) fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    let Self { partition, gtr, .. } = self;
    let MarginalUpdate {
      node_states,
      edges,
      log_lh,
    } = partition.marginal_update(&gtr, graph, branch_lengths, &())?;
    Ok((
      Self {
        partition,
        gtr,
        node_states,
        edges,
      },
      log_lh,
    ))
  }
}

#[derive(Debug, Serialize)]
pub struct AncestralOutput {
  #[serde(skip)]
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub mask: Vec<bool>,
  #[serde(skip)]
  pub emitted_nodes: Vec<GraphNodeKey>,
}
