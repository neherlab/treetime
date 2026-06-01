use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::ClockParams;
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::partition::timetree::GraphTimetree;
use crate::partition::traits::PartitionTimetreeAll;
use crate::payload::timetree::EdgeTimetree;
use crate::payload::timetree::NodeTimetree;
use crate::timetree::refinement::{self, RefinementParams};
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_distribution::Distribution;

pub fn run_refinement_iteration(
  args: &TreetimeTimetreeArgs,
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  clock_model: &mut ClockModel,
  clock_params: &ClockParams,
  branch_params: &BranchPointOptimizationParams,
  coalescent_tc: Option<&Distribution>,
) -> Result<(usize, usize), Report> {
  let params = RefinementParams {
    relax: args.relax.clone(),
    resolve_polytomies: args.resolve_polytomies,
    clock_rate: args.clock_rate,
    no_indels: args.no_indels,
  };
  refinement::run_refinement_iteration(&params, graph, partitions, clock_model, clock_params, branch_params, coalescent_tc)
}
