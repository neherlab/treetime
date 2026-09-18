//! Shared output-planning and tree-ordering steps for the server's command orchestrations.
//!
//! The server always requests a command's default `--output-all` file set into one directory and
//! applies the default topology ordering before writing, matching the command-line client's default
//! behavior. These helpers build the fixed output plan and topology spec so each command reads them
//! the same way.

use app_output::output_plan::{CommandKind, OutputPlanRequest, OutputSelection, ResolvedOutputs, plan};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::clock::find_best_root::params::{RerootMethod, RerootSpec};
use treetime_graph::topology_order::{TopologyOrderPreset, TopologyOrderSpec, TopologyOrderTargetAggregate};

/// Resolve a command's `--output-all` file set into concrete paths under `outdir`.
///
/// The server requests the command's default set into one directory and, for the few requests that
/// carry a per-file destination (mugration confidence CSV, timetree tracelog), passes it as a
/// non-tree override. It never sets a restricting selection or a non-default annotation style.
pub fn output_plan(
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

/// Resolve a command's default `--output-all` file set with no per-file overrides.
pub fn default_output_plan(command: CommandKind, outdir: &Path) -> Result<ResolvedOutputs, Report> {
  output_plan(command, outdir, BTreeMap::new())
}

/// Topology ordering applied before writing: descendant-count ladderization, the command-line
/// client's default when neither `--ladderize` nor `--topology-order` is given.
pub fn default_topology_order() -> TopologyOrderSpec {
  TopologyOrderSpec {
    preset: TopologyOrderPreset::DescendantCount,
    target_order: vec![],
    target_aggregate: TopologyOrderTargetAggregate::Mean,
  }
}

/// Requested reroot policy from a method plus a tip list: a non-empty tip list reroots on the MRCA of
/// the tips, otherwise the method applies (least-squares by default).
pub fn reroot_spec(reroot: Option<RerootMethod>, reroot_tips: &[String]) -> RerootSpec {
  if reroot_tips.is_empty() {
    RerootSpec::Method(reroot.unwrap_or_default())
  } else {
    RerootSpec::Tips(reroot_tips.to_vec())
  }
}
