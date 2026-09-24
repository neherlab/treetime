use app_output::output_plan::{CommandKind, OutputPlanRequest, OutputSelection, ResolvedOutputs, plan};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
use treetime_graph::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};

pub(crate) fn default_output_plan(command: CommandKind, outdir: &Path) -> Result<ResolvedOutputs, Report> {
  output_plan(command, outdir, BTreeMap::new())
}

pub(crate) fn output_plan(
  command: CommandKind,
  outdir: &Path,
  non_tree_overrides: BTreeMap<OutputSelection, PathBuf>,
) -> Result<ResolvedOutputs, Report> {
  plan(&OutputPlanRequest {
    command,
    output_all: Some(outdir.to_path_buf()),
    nwk_styles: vec![],
    selection: vec![],
    tree_overrides: BTreeMap::new(),
    non_tree_overrides,
  })
}

pub(crate) fn default_topology_order() -> TopologyOrderSpec {
  TopologyOrderSpec {
    preset: TopologyOrderPreset::DescendantCount,
    target_order: vec![],
    target_aggregate: TopologyOrderTargetAggregate::Mean,
  }
}

pub(crate) fn reroot_spec(reroot: Option<RerootMethod>, reroot_tips: &[String]) -> RerootSpec {
  if reroot_tips.is_empty() {
    RerootSpec::Method(reroot.unwrap_or_default())
  } else {
    RerootSpec::Tips(reroot_tips.to_vec())
  }
}
