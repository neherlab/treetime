use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use parking_lot::RwLock;
use serde::Serialize;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;

#[derive(Serialize)]
pub struct OptimizeGraphData {
  pub gtr: GTR,
  pub model_name: GtrModelName,
  pub sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  pub dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>>,
}

impl OptimizeGraphData {
  pub fn new(
    gtr: GTR,
    model_name: GtrModelName,
    sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
    dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>>,
  ) -> Self {
    Self {
      gtr,
      model_name,
      sparse_partitions,
      dense_partitions,
    }
  }
}

/// Per-node optimize output as a value: the name and input branch support the output writers read.
#[derive(Debug, Clone, Serialize)]
pub struct OptimizeNodeOut {
  pub name: Option<String>,
  pub confidence: Option<f64>,
}

/// Per-edge optimize output as a value: the optimized branch length the output writers read.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct EdgeOut {
  pub branch_length: Option<f64>,
}

/// Branch-length optimization result as a value.
///
/// `edges` holds the optimized per-edge branch lengths keyed by edge id, and the sequence
/// partitions plus the substitution model sit alongside. `graph` carries the tree the output
/// writers still read from; the writers move onto the result value in a later step, after which
/// `graph` and the partition-in-graph go away.
#[derive(Serialize)]
pub struct OptimizeResult {
  #[serde(skip)]
  pub graph: GraphAncestral<OptimizeGraphData>,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, OptimizeNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub gtr: GTR,
  #[serde(skip)]
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub sparse_partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  #[serde(skip)]
  pub dense_partitions: Vec<Arc<RwLock<PartitionMarginalDense>>>,
}
