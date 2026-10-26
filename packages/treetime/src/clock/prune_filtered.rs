use crate::progress::{LogSink, RunWarning, RunWarningKind};
use crate::progress_info;
use crate::reroot::params::RerootSpec;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeSet;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::remove_leaves::remove_leaves;
use treetime_utils::fmt::string::{counted, name_list};
use treetime_utils::make_error;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PruneFilteredTips {
  pub keep: BTreeSet<GraphNodeKey>,
}

pub(crate) struct FilteredTips {
  pub removed: BTreeSet<GraphNodeKey>,
  pub kept: BTreeSet<GraphNodeKey>,
}

impl PruneFilteredTips {
  pub(crate) fn split(&self, outliers: BTreeSet<GraphNodeKey>) -> FilteredTips {
    let (kept, removed) = outliers.into_iter().partition(|tip| self.keep.contains(tip));
    FilteredTips { removed, kept }
  }
}

pub(crate) fn remove_filtered_tips(
  graph: &mut Graph,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  tips: &FilteredTips,
  reroot_spec: Option<&RerootSpec>,
  names: &NodeNames,
  log: &dyn LogSink,
) -> Result<(), Report> {
  if let Some(RerootSpec::Tips(outgroup)) = reroot_spec
    && let Some(&tip) = outgroup.iter().find(|tip| tips.removed.contains(tip))
  {
    return make_error!(
      "The clock filter flagged the outgroup tip {}, and --prune-filtered-tips would remove it. Name it in --keep-tips to keep it in the tree",
      &names[tip]
    );
  }
  if !tips.kept.is_empty() {
    progress_info!(
      log,
      "Keeping {} that the clock filter flagged and the keep list names",
      counted(tips.kept.len(), "tip")
    );
  }
  if tips.removed.is_empty() {
    return Ok(());
  }
  remove_leaves(graph, branch_lengths, &tips.removed)?;
  let removed = tips
    .removed
    .iter()
    .map(|&tip| names[tip].to_owned())
    .sorted()
    .collect_vec();
  log.warning(&RunWarning {
    kind: RunWarningKind::FilteredTipsRemoved,
    message: format!(
      "Removed {} that the clock filter flagged from the tree ({})",
      counted(removed.len(), "tip"),
      name_list(&removed)
    ),
    names: removed,
  });
  Ok(())
}
