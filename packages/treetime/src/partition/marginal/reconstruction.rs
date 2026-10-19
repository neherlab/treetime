use crate::alphabet::alphabet::Alphabet;
use crate::error::OperationError;
use crate::gtr::gtr::GTR;
use crate::gtr::refinement::{RateRangeEdge, optimize_gtr_rate};
use crate::partition::marginal::dense::partition::{DenseMarginalEdges, PartitionMarginalDense};
use crate::partition::marginal::sample::SampleMode;
use crate::partition::marginal::shared::reconcile::{live_node_keys, reconcile_node_states};
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::partition::marginal::sparse::mutations::sparse_edge_mutations;
use crate::partition::marginal::sparse::partition::{PartitionMarginalSparse, SparseMarginalEdges};
use crate::partition::marginal::sparse::reroot::reroot_sparse;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::DenseNodeState;
use crate::partition::storage::sparse::SparseNodeState;
use crate::progress::LogSink;
use crate::seq::indel::InDel;
use crate::seq::mutation::{
  Mutation, MutationTrack, SequenceMutations, Sub, combine_edge_mutations, emit_sequences, stream_sequence_mutations,
};
use crate::seq::overlay::SeqOverlay;
use crate::seq::sink::{SeqSink, SeqTrack};
use eyre::Report;
use rand::Rng;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootResult;
use treetime_primitives::{LogLh, Seq};

#[derive(Clone, Debug)]
pub enum MarginalReconstruction {
  Dense(DenseReconstruction),
  Sparse(SparseReconstruction),
}

macro_rules! each_kind {
  ($self:expr, |$reconstruction:ident| $body:expr) => {
    match $self {
      Self::Dense($reconstruction) => $body,
      Self::Sparse($reconstruction) => $body,
    }
  };
}

impl MarginalReconstruction {
  pub(crate) fn view(&self) -> ReconstructionRef<'_> {
    match self {
      Self::Dense(reconstruction) => ReconstructionRef::Dense(reconstruction),
      Self::Sparse(reconstruction) => ReconstructionRef::Sparse(reconstruction),
    }
  }

  pub(crate) fn gtr(&self) -> &GTR {
    self.view().gtr()
  }

  pub fn alphabet(&self) -> &Alphabet {
    self.view().alphabet()
  }

  pub(crate) fn graph_log_lh(&self, graph: &Graph) -> Result<LogLh, Report> {
    self.view().graph_log_lh(graph)
  }

  pub(crate) fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    Ok(match self {
      Self::Dense(reconstruction) => {
        let (reconstruction, log_lh) = reconstruction.marginal_update(graph, branch_lengths)?;
        (Self::Dense(reconstruction), log_lh)
      },
      Self::Sparse(reconstruction) => {
        let (reconstruction, log_lh) = reconstruction.marginal_update(graph, branch_lengths)?;
        (Self::Sparse(reconstruction), log_lh)
      },
    })
  }

  pub(crate) fn sample_sequences(
    &self,
    graph: &Graph,
    sample_mode: SampleMode,
    rng: &mut dyn Rng,
  ) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
    each_kind!(self, |reconstruction| reconstruction.partition.sample_sequences(
      graph,
      &reconstruction.node_states,
      sample_mode,
      rng
    ))
  }

  pub(crate) fn node_sequence(&self, graph: &Graph, impute: bool, node_key: GraphNodeKey) -> Result<Seq, Report> {
    self.view().node_sequence(graph, impute, node_key)
  }

  pub(crate) fn stream_sequences(
    &self,
    graph: &Graph,
    impute: bool,
    track: SeqTrack<'_>,
    include_leaves: bool,
    report_unknown: bool,
    sink: Option<&mut (dyn SeqSink + '_)>,
  ) -> Result<SequenceMutations, OperationError> {
    self
      .view()
      .stream_sequences(graph, impute, track, include_leaves, report_unknown, sink)
  }

  pub(crate) fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    self.view().edge_indels(edge_key)
  }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum ReconstructionRef<'a> {
  Dense(&'a DenseReconstruction),
  Sparse(&'a SparseReconstruction),
}

impl<'a> ReconstructionRef<'a> {
  pub(crate) fn gtr(self) -> &'a GTR {
    each_kind!(self, |reconstruction| &reconstruction.gtr)
  }

  pub(crate) fn length(self) -> usize {
    each_kind!(self, |reconstruction| reconstruction.partition.length)
  }

  pub(crate) fn alphabet(self) -> &'a Alphabet {
    each_kind!(self, |reconstruction| &reconstruction.partition.alphabet)
  }

  pub(crate) fn graph_log_lh(self, graph: &Graph) -> Result<LogLh, Report> {
    let root_key = graph.root_key()?;
    Ok(each_kind!(self, |reconstruction| reconstruction
      .partition
      .get_log_lh(&reconstruction.node_states, root_key)))
  }

  pub(crate) fn node_sequence(self, graph: &Graph, impute: bool, node_key: GraphNodeKey) -> Result<Seq, Report> {
    match self {
      Self::Dense(reconstruction) => {
        reconstruction
          .partition
          .node_sequence(graph, &reconstruction.node_states, impute, node_key)
      },
      Self::Sparse(reconstruction) => reconstruction.partition.node_sequence(
        graph,
        &reconstruction.node_states,
        &reconstruction.edges.forward,
        impute,
        node_key,
      ),
    }
  }

  pub(crate) fn stream_sequences(
    self,
    graph: &Graph,
    impute: bool,
    track: SeqTrack<'_>,
    include_leaves: bool,
    report_unknown: bool,
    sink: Option<&mut (dyn SeqSink + '_)>,
  ) -> Result<SequenceMutations, OperationError> {
    let node_sequence = |node_key| self.node_sequence(graph, impute, node_key);
    let Self::Sparse(reconstruction) = self else {
      return stream_sequence_mutations(
        graph,
        self.alphabet(),
        track,
        include_leaves,
        report_unknown,
        node_sequence,
        |edge_key| self.edge_indels(edge_key),
        sink,
      );
    };
    if let Some(sink) = sink {
      emit_sequences(graph, track, include_leaves, node_sequence, sink)?;
    }
    let root_key = graph.root_key().map_err(OperationError::classify)?;
    let root_sequence = node_sequence(root_key).map_err(OperationError::classify)?;
    let edge_mutations = sparse_edge_mutations(
      &reconstruction.partition,
      graph,
      &reconstruction.node_states,
      &reconstruction.edges.forward,
      impute,
      report_unknown,
      &track.mutation_track(),
    )
    .map_err(OperationError::classify)?;
    Ok(SequenceMutations {
      root_sequence,
      edge_mutations,
    })
  }

  pub(crate) fn extract_ancestral_sequence(self, node_key: GraphNodeKey) -> SeqOverlay {
    each_kind!(self, |reconstruction| reconstruction
      .partition
      .extract_ancestral_sequence(&reconstruction.node_states, node_key))
  }

  pub(crate) fn root_sequence(self, graph: &Graph) -> Result<Seq, Report> {
    each_kind!(self, |reconstruction| reconstruction
      .partition
      .root_sequence(&reconstruction.node_states, graph))
  }

  pub(crate) fn edge_subs(self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .edge_subs(&reconstruction.node_states, graph, edge_key),
      Self::Sparse(reconstruction) => reconstruction
        .partition
        .edge_subs(&reconstruction.edges.estimates, edge_key),
    }
  }

  pub(crate) fn edge_sub_count(self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Option<usize>, Report> {
    match self {
      Self::Dense(_) => Ok(Some(self.edge_subs(graph, edge_key)?.len())),
      Self::Sparse(reconstruction) => Ok(reconstruction.edges.estimates.get(&edge_key).map(Vec::len)),
    }
  }

  pub(crate) fn edge_indels(self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .edge_indels(&reconstruction.edges.estimates, edge_key),
      Self::Sparse(reconstruction) => reconstruction.partition.edge_indels(edge_key),
    }
  }

  pub(crate) fn edge_mutations(
    self,
    graph: &Graph,
    edge_key: GraphEdgeKey,
    track: &MutationTrack,
  ) -> Result<Vec<Mutation>, Report> {
    combine_edge_mutations(self.edge_subs(graph, edge_key)?, &self.edge_indels(edge_key), track)
  }

  pub(crate) fn edge_effective_length(self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    match self {
      Self::Dense(reconstruction) => {
        reconstruction
          .partition
          .edge_effective_length(&reconstruction.node_states, graph, edge_key)
      },
      Self::Sparse(reconstruction) => reconstruction.partition.edge_effective_length(graph, edge_key),
    }
  }

  pub(crate) fn create_edge_contribution(self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    match self {
      Self::Dense(reconstruction) => reconstruction.partition.create_edge_contribution(
        &reconstruction.gtr,
        &reconstruction.edges.backward,
        &reconstruction.edges.forward,
        edge_key,
      ),
      Self::Sparse(reconstruction) => reconstruction.partition.create_edge_contribution(
        &reconstruction.gtr,
        &reconstruction.edges.backward,
        &reconstruction.edges.forward,
        edge_key,
      ),
    }
  }

  pub(crate) fn edge_has_fitch_mutations(self, edge_key: GraphEdgeKey) -> bool {
    match self {
      Self::Dense(_) => false,
      Self::Sparse(reconstruction) => reconstruction
        .partition
        .obs_edges
        .get(&edge_key)
        .is_some_and(|edge| !edge.fitch_subs().is_empty() || !edge.indels.is_empty()),
    }
  }

  pub(crate) fn edge_indel_count(self, edge_key: GraphEdgeKey) -> usize {
    match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .edge_indel_count(&reconstruction.edges.estimates, edge_key),
      Self::Sparse(reconstruction) => reconstruction.partition.edge_indel_count(edge_key),
    }
  }
}

#[derive(Clone, Debug)]
pub struct SparseReconstruction {
  pub(crate) partition: PartitionMarginalSparse,
  pub(crate) gtr: GTR,
  pub(crate) node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  pub(crate) edges: SparseMarginalEdges,
}

impl SparseReconstruction {
  pub(crate) fn seeded(
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

  pub(crate) fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    let Self {
      partition,
      gtr,
      node_states,
      edges,
    } = self;
    drop(edges);
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

  #[must_use]
  pub(crate) fn with_total_length(self, total_length: usize) -> Self {
    Self {
      partition: self.partition.with_total_length(total_length),
      ..self
    }
  }

  pub(crate) fn with_estimated_rate(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    log: &dyn LogSink,
  ) -> Result<Self, Report> {
    let (gtr, _, _) = optimize_gtr_rate(
      &self.partition,
      self.gtr,
      &self.node_states,
      graph,
      branch_lengths,
      RateRangeEdge::LowerCostEnd,
      log,
    )?;
    Ok(Self::seeded(self.partition, gtr, self.node_states))
  }

  pub(crate) fn apply_reroot(self, changes: &RerootResult) -> Result<Self, Report> {
    reroot_sparse(self.partition, self.gtr, self.node_states, changes)
  }

  #[must_use]
  pub(crate) fn reconcile_topology(self, graph: &Graph) -> Self {
    let live_nodes = live_node_keys(graph);
    let mut partition = self.partition;
    partition.reconcile_topology(graph);
    Self::seeded(
      partition,
      self.gtr,
      reconcile_node_states(self.node_states, &live_nodes, SparseNodeState::empty),
    )
  }

  #[must_use]
  pub(crate) fn reconcile_node_states(self, graph: &Graph) -> Self {
    Self::seeded(
      self.partition,
      self.gtr,
      reconcile_node_states(self.node_states, &live_node_keys(graph), SparseNodeState::empty),
    )
  }
}

#[derive(Clone, Debug)]
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

  #[must_use]
  pub(crate) fn with_total_length(self, total_length: usize) -> Self {
    Self {
      partition: self.partition.with_total_length(total_length),
      ..self
    }
  }

  #[must_use]
  pub(crate) fn reseeded(self) -> Self {
    Self::seeded(self.partition, self.gtr)
  }

  pub(crate) fn with_estimated_rate(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    log: &dyn LogSink,
  ) -> Result<Self, Report> {
    let (gtr, _, _) = optimize_gtr_rate(
      &self.partition,
      self.gtr,
      &self.node_states,
      graph,
      branch_lengths,
      RateRangeEdge::LowerCostEnd,
      log,
    )?;
    Ok(Self::seeded(self.partition, gtr))
  }

  pub(crate) fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<(Self, LogLh), Report> {
    let Self {
      partition,
      gtr,
      node_states,
      edges,
    } = self;
    drop(node_states);
    drop(edges);
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
