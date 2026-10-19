use crate::branch_lengths::branch_lengths_or_zero;
use crate::partition::set::PartitionSet;
use crate::progress::LogSink;
use crate::progress_info;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::reroot::{EdgeMergeInfo, RerootResult};

pub(crate) enum BranchModel {
  Input,
  Marginal(PartitionSet),
}

impl BranchModel {
  pub(crate) fn is_marginal(&self) -> bool {
    matches!(self, Self::Marginal(_))
  }

  pub(crate) fn marginal_update(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  ) -> Result<Self, Report> {
    Ok(match self {
      Self::Input => Self::Input,
      Self::Marginal(partitions) => Self::Marginal(partitions.marginal_update(graph, branch_lengths)?.0),
    })
  }

  pub(crate) fn reconcile_topology(self, graph: &Graph, merges: &[EdgeMergeInfo]) -> Result<Self, Report> {
    Ok(match self {
      Self::Input => Self::Input,
      Self::Marginal(partitions) => Self::Marginal(partitions.reconcile_topology(graph, merges)?),
    })
  }

  pub(crate) fn apply_reroot(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    reroot: &RerootResult,
    log: &dyn LogSink,
  ) -> Result<Self, Report> {
    let Self::Marginal(partitions) = self else {
      return Ok(self);
    };
    progress_info!(log, "Applying reroot changes to the marginal reconstruction");
    let partitions = partitions
      .apply_reroot(reroot)
      .wrap_err("Failed to apply reroot changes to the partitions")?
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))
      .wrap_err("Failed to update marginal after reroot")?
      .0;
    Ok(Self::Marginal(partitions))
  }
}
