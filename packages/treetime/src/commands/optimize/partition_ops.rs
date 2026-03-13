use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::seq::mutation::Sub;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;

pub trait PartitionOptimizeOps: Send + Sync {
  fn sequence_length(&self) -> usize;

  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;

  fn edge_subs(&self, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;
}

pub type PartitionOptimizeVec = Vec<Arc<RwLock<dyn PartitionOptimizeOps>>>;
