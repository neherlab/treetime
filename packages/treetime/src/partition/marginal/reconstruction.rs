use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::{DenseMarginalEdges, PartitionMarginalDense};
use crate::partition::marginal::sample::SampleMode;
use crate::partition::marginal::sequences::TipStates;
use crate::partition::marginal::sequences::{ReconstructedSequences, emitted_nodes};
use crate::partition::marginal::shared::reconcile::{live_node_keys, reconcile_node_states};
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::partition::marginal::sparse::partition::{PartitionMarginalSparse, SparseMarginalEdges};
use crate::partition::marginal::sparse::reroot::reroot_sparse;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::DenseNodeState;
use crate::partition::storage::sparse::SparseNodeState;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
use eyre::Report;
use rand::RngCore;
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootResult;
use treetime_primitives::{AsciiChar, LogLh, Seq};

#[derive(Clone, Debug, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum MarginalReconstruction {
  Dense(DenseReconstruction),
  Sparse(SparseReconstruction),
}

impl MarginalReconstruction {
  pub(crate) fn gtr(&self) -> &GTR {
    match self {
      Self::Dense(reconstruction) => &reconstruction.gtr,
      Self::Sparse(reconstruction) => &reconstruction.gtr,
    }
  }

  pub(crate) fn gtr_mut(&mut self) -> &mut GTR {
    match self {
      Self::Dense(reconstruction) => &mut reconstruction.gtr,
      Self::Sparse(reconstruction) => &mut reconstruction.gtr,
    }
  }

  pub(crate) fn sequence_length(&self) -> usize {
    match self {
      Self::Dense(reconstruction) => reconstruction.partition.length,
      Self::Sparse(reconstruction) => reconstruction.partition.length,
    }
  }

  pub(crate) fn ambiguous_char(&self) -> AsciiChar {
    match self {
      Self::Dense(reconstruction) => reconstruction.partition.alphabet.unknown(),
      Self::Sparse(reconstruction) => reconstruction.partition.alphabet.unknown(),
    }
  }

  pub(crate) fn graph_log_lh(&self, graph: &Graph) -> Result<LogLh, Report> {
    let root_key = graph.root_key()?;
    Ok(match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .get_log_lh(&reconstruction.node_states, root_key),
      Self::Sparse(reconstruction) => reconstruction
        .partition
        .get_log_lh(&reconstruction.node_states, root_key),
    })
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

  pub(crate) fn reconstruct_sequences(
    &self,
    graph: &Graph,
    tips: TipStates,
    sample_mode: SampleMode,
    rng: &mut dyn RngCore,
  ) -> Result<ReconstructedSequences, Report> {
    let sampled = match self {
      Self::Dense(reconstruction) => {
        reconstruction
          .partition
          .sample_sequences(graph, &reconstruction.node_states, sample_mode, rng)?
      },
      Self::Sparse(reconstruction) => {
        reconstruction
          .partition
          .sample_sequences(graph, &reconstruction.node_states, sample_mode, rng)?
      },
    };
    Ok(ReconstructedSequences {
      sampled,
      emitted_nodes: emitted_nodes(graph, tips.include_leaves)?,
    })
  }

  pub(crate) fn node_sequence(&self, graph: &Graph, impute: bool, node_key: GraphNodeKey) -> Result<Seq, Report> {
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

  pub(crate) fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .extract_ancestral_sequence(&reconstruction.node_states, node_key),
      Self::Sparse(reconstruction) => reconstruction
        .partition
        .extract_ancestral_sequence(&reconstruction.node_states, node_key),
    }
  }

  pub fn root_sequence(&self, graph: &Graph) -> Result<Seq, Report> {
    match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .root_sequence(&reconstruction.node_states, graph),
      Self::Sparse(reconstruction) => reconstruction
        .partition
        .root_sequence(&reconstruction.node_states, graph),
    }
  }

  pub fn edge_subs(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report> {
    match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .edge_subs(&reconstruction.node_states, graph, edge_key),
      Self::Sparse(reconstruction) => reconstruction
        .partition
        .edge_subs(&reconstruction.edges.estimates, edge_key),
    }
  }

  pub(crate) fn edge_sub_count(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<Option<usize>, Report> {
    match self {
      Self::Dense(_) => Ok(Some(self.edge_subs(graph, edge_key)?.len())),
      Self::Sparse(reconstruction) => Ok(reconstruction.edges.estimates.get(&edge_key).map(Vec::len)),
    }
  }

  pub(crate) fn edge_indels(&self, edge_key: GraphEdgeKey) -> Vec<InDel> {
    match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .edge_indels(&reconstruction.edges.estimates, edge_key),
      Self::Sparse(reconstruction) => reconstruction.partition.edge_indels(edge_key),
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

  pub(crate) fn edge_effective_length(&self, graph: &Graph, edge_key: GraphEdgeKey) -> Result<usize, Report> {
    match self {
      Self::Dense(reconstruction) => {
        reconstruction
          .partition
          .edge_effective_length(&reconstruction.node_states, graph, edge_key)
      },
      Self::Sparse(reconstruction) => reconstruction.partition.edge_effective_length(graph, edge_key),
    }
  }

  pub(crate) fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    match self {
      Self::Dense(reconstruction) => Ok(reconstruction.partition.create_edge_contribution(
        &reconstruction.gtr,
        &reconstruction.edges.backward,
        &reconstruction.edges.forward,
        edge_key,
      )),
      Self::Sparse(reconstruction) => reconstruction.partition.create_edge_contribution(
        &reconstruction.gtr,
        &reconstruction.edges.backward,
        &reconstruction.edges.forward,
        edge_key,
      ),
    }
  }

  pub(crate) fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    match self {
      Self::Dense(reconstruction) => reconstruction
        .partition
        .edge_indel_count(&reconstruction.edges.estimates, edge_key),
      Self::Sparse(reconstruction) => reconstruction.partition.edge_indel_count(edge_key),
    }
  }

  pub(crate) fn apply_reroot(self, changes: &RerootResult) -> Result<Self, Report> {
    Ok(match self {
      Self::Dense(reconstruction) => Self::Dense(DenseReconstruction::seeded(
        reconstruction.partition,
        reconstruction.gtr,
      )),
      Self::Sparse(reconstruction) => Self::Sparse(reroot_sparse(
        reconstruction.partition,
        reconstruction.gtr,
        reconstruction.node_states,
        changes,
      )?),
    })
  }

  #[must_use]
  pub(crate) fn reconcile_topology(self, graph: &Graph) -> Self {
    match self {
      Self::Dense(reconstruction) => Self::Dense(DenseReconstruction::seeded(
        reconstruction.partition,
        reconstruction.gtr,
      )),
      Self::Sparse(reconstruction) => {
        let live_nodes = live_node_keys(graph);
        let mut partition = reconstruction.partition;
        partition.reconcile_topology(graph);
        Self::Sparse(SparseReconstruction::seeded(
          partition,
          reconstruction.gtr,
          reconcile_node_states(reconstruction.node_states, &live_nodes, SparseNodeState::empty),
        ))
      },
    }
  }
}

#[derive(Clone, Debug, Serialize)]
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
