use crate::partition::marginal::dense::reroot::reroot_dense;
use crate::partition::marginal::sparse::reroot::reroot_sparse;
use crate::partition::optimize::contribution::OptimizationContribution;
use crate::partition::storage::dense::DenseNodeState;
use crate::partition::storage::sparse::SparseNodeState;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::partition::traits::PartitionOptimizeOps;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootChanges;

impl PartitionOptimizeOps for PartitionTimetree {
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report> {
    match self {
      Self::Dense(family) => family.readout().create_edge_contribution(edge_key),
      Self::Sparse(family) => family.readout().create_edge_contribution(edge_key),
    }
  }

  fn edge_indel_count(&self, edge_key: GraphEdgeKey) -> usize {
    match self {
      Self::Dense(family) => family.readout().edge_indel_count(edge_key),
      Self::Sparse(family) => family.readout().edge_indel_count(edge_key),
    }
  }
}

impl PartitionTimetree {
  /// Apply a reroot to this partition's observations and result maps at the structural operation.
  pub fn apply_reroot(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    match self {
      Self::Dense(family) => reroot_dense(family, changes),
      Self::Sparse(family) => reroot_sparse(family, changes),
    }
  }

  /// Ensure the partition observations and node states have entries for all nodes and edges in the
  /// graph after a topology change (polytomy resolution), and drop entries no longer in the graph.
  /// The subsequent marginal update pass recomputes values; the stale messages and estimates are
  /// dropped here.
  pub fn reconcile_topology(&mut self, graph: &Graph) {
    match self {
      Self::Dense(family) => {
        let node_keys: Vec<GraphNodeKey> = graph.get_nodes().iter().map(|node| node.read_arc().key()).collect();
        for &key in &node_keys {
          family.node_states.entry(key).or_insert_with(DenseNodeState::empty);
        }
        family.node_states.retain(|key, _| node_keys.contains(key));
        family.backward.clear();
        family.forward.clear();
        family.estimates.clear();
      },
      Self::Sparse(family) => {
        family.partition.reconcile_topology(graph);
        let node_keys: Vec<GraphNodeKey> = graph.get_nodes().iter().map(|node| node.read_arc().key()).collect();
        for &key in &node_keys {
          family.node_states.entry(key).or_insert_with(SparseNodeState::empty);
        }
        family.node_states.retain(|key, _| node_keys.contains(key));
        family.backward.clear();
        family.forward.clear();
        family.estimates.clear();
      },
    }
  }
}
