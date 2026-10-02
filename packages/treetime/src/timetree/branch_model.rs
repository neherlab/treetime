use crate::ancestral::marginal::branch_lengths_or_zero;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::progress::ProgressSink;
use crate::progress_info;
use eyre::{Report, WrapErr};
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::reroot::{RerootChanges, RerootResult};

#[expect(
  clippy::large_enum_variant,
  reason = "one value per run moves between steps by value; it is never stored in a collection"
)]
pub(crate) enum BranchModel {
  Input,
  Marginal(PartitionTimetree),
}

impl BranchModel {
  pub(crate) fn sequence_length(&self) -> usize {
    match self {
      Self::Input => 0,
      Self::Marginal(partition) => partition.sequence_length(),
    }
  }

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
      Self::Marginal(partition) => Self::Marginal(partition.marginal_update(graph, branch_lengths)?),
    })
  }

  #[must_use]
  pub(crate) fn reconcile_topology(self, graph: &Graph) -> Self {
    match self {
      Self::Input => Self::Input,
      Self::Marginal(partition) => Self::Marginal(partition.reconcile_topology(graph)),
    }
  }

  pub(crate) fn apply_reroot(
    self,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    reroot: &RerootResult,
    progress: &dyn ProgressSink,
  ) -> Result<Self, Report> {
    let Self::Marginal(partition) = self else {
      return Ok(self);
    };
    let changes = RerootChanges {
      edge_split: reroot.edge_split.clone(),
      edge_merge: reroot.edge_merge.clone(),
      inverted_edge_keys: reroot.inverted_edge_keys.clone(),
    };
    progress_info!(progress, "Applying reroot changes to 1 partitions");
    let partition = partition
      .apply_reroot(&changes)
      .wrap_err("Failed to apply reroot changes to partition")?
      .marginal_update(graph, &branch_lengths_or_zero(branch_lengths))
      .wrap_err("Failed to update marginal after reroot")?;
    Ok(Self::Marginal(partition))
  }
}
