use crate::progress::{LogSink, RunWarning, RunWarningKind};
use crate::progress_info;
use crate::reroot::params::RerootSpec;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeSet;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::remove_leaves::remove_leaves;
use treetime_primitives::date::DateConstraint;
use treetime_utils::fmt::string::{counted, name_list};
use treetime_utils::make_error;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PruneTips {
  pub filtered: bool,
  pub undated: bool,
  pub keep: BTreeSet<GraphNodeKey>,
}

pub(crate) struct TipRemoval {
  pub flagged: BTreeSet<GraphNodeKey>,
  pub undated: BTreeSet<GraphNodeKey>,
  pub kept_flagged: BTreeSet<GraphNodeKey>,
  pub kept_undated: BTreeSet<GraphNodeKey>,
}

impl PruneTips {
  pub fn is_active(&self) -> bool {
    self.filtered || self.undated
  }

  pub(crate) fn split(
    &self,
    flagged: BTreeSet<GraphNodeKey>,
    graph: &Graph,
    dates: &NodeMap<DateConstraint>,
  ) -> TipRemoval {
    let flagged = if self.filtered { flagged } else { BTreeSet::new() };
    let undated = if self.undated {
      undated_tips(graph, dates)
    } else {
      BTreeSet::new()
    };
    let (kept_flagged, flagged) = flagged.into_iter().partition(|tip| self.keep.contains(tip));
    let (kept_undated, undated) = undated.into_iter().partition(|tip| self.keep.contains(tip));
    TipRemoval {
      flagged,
      undated,
      kept_flagged,
      kept_undated,
    }
  }
}

pub(crate) fn remove_tips(
  graph: &mut Graph,
  branch_lengths: &mut EdgeMap<Option<f64>>,
  tips: &TipRemoval,
  reroot_spec: Option<&RerootSpec>,
  names: &NodeNames,
  log: &dyn LogSink,
) -> Result<(), Report> {
  if let Some(RerootSpec::Tips(outgroup)) = reroot_spec {
    if let Some(&tip) = outgroup.iter().find(|tip| tips.flagged.contains(tip)) {
      return make_error!(
        "The clock filter flagged the outgroup tip {}, and --prune-filtered-tips would remove it. Name it in --keep-tips to keep it in the tree",
        &names[tip]
      );
    }
    if let Some(&tip) = outgroup.iter().find(|tip| tips.undated.contains(tip)) {
      return make_error!(
        "The outgroup tip {} has no usable date, and --prune-undated-tips would remove it. Name it in --keep-tips to keep it in the tree",
        &names[tip]
      );
    }
  }
  if !tips.kept_flagged.is_empty() {
    progress_info!(
      log,
      "Keeping {} that the clock filter flagged and the keep list names",
      counted(tips.kept_flagged.len(), "tip")
    );
  }
  if !tips.kept_undated.is_empty() {
    progress_info!(
      log,
      "Keeping {} without a usable date that the keep list names",
      counted(tips.kept_undated.len(), "tip")
    );
  }
  let removed: BTreeSet<GraphNodeKey> = tips.flagged.union(&tips.undated).copied().collect();
  if removed.is_empty() {
    return Ok(());
  }
  remove_leaves(graph, branch_lengths, &removed)?;
  warn_removed(
    &tips.flagged,
    RunWarningKind::FilteredTipsRemoved,
    "that the clock filter flagged",
    names,
    log,
  );
  warn_removed(
    &tips.undated,
    RunWarningKind::UndatedTipsRemoved,
    "without a usable date",
    names,
    log,
  );
  Ok(())
}

fn undated_tips(graph: &Graph, dates: &NodeMap<DateConstraint>) -> BTreeSet<GraphNodeKey> {
  graph
    .get_leaves()
    .map(|leaf| leaf.key())
    .filter(|&key| !dates.contains_key(key))
    .collect()
}

fn warn_removed(
  tips: &BTreeSet<GraphNodeKey>,
  kind: RunWarningKind,
  description: &str,
  names: &NodeNames,
  log: &dyn LogSink,
) {
  if tips.is_empty() {
    return;
  }
  let removed = tips.iter().map(|&tip| names[tip].to_owned()).sorted().collect_vec();
  log.warning(&RunWarning {
    kind,
    message: format!(
      "Removed {} {description} from the tree ({})",
      counted(removed.len(), "tip"),
      name_list(&removed)
    ),
    names: removed,
  });
}
