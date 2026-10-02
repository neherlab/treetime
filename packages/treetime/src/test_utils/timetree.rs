use crate::alphabet::alphabet::AlphabetName;
use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::RerootSpec;
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
use crate::gtr::get_gtr::GtrModelName;
use crate::optimize::params::BranchLengthMode;
use crate::seq::gap_fill::GapFill;
use crate::timetree::inference::bad_branches::{derive_bad_branches, undated_leaves};
use crate::timetree::inference::time_inference::{
  BranchLikelihood, NodePosterior, TimeBackward, TimeInference, likely_times,
};
use crate::timetree::params::TimeMarginalMode;
use crate::timetree::pipeline::TimetreeParams;
use eyre::Report;
use treetime_graph::graph::Graph;

pub(crate) fn constraint_coalescent_node_times(
  graph: &Graph,
  constraints: &DateConstraints,
) -> Result<CoalescentNodeTimes, Report> {
  let bad_branches = derive_bad_branches(graph, constraints, &undated_leaves(graph, constraints))?;
  let times = likely_times(graph, constraints, None)?
    .into_iter()
    .map(|(key, time_dist_likely)| {
      let entry = CoalescentNodeTime {
        time: None,
        time_dist_likely,
        bad_branch: bad_branches[&key],
      };
      (key, entry)
    })
    .collect();
  Ok(times)
}

pub(crate) fn empty_time_inference(graph: &Graph) -> TimeInference {
  TimeInference {
    bad_branches: graph.get_nodes().map(|node| (node.key(), false)).collect(),
    branches: graph
      .get_edges()
      .map(|edge| {
        let branch = BranchLikelihood {
          distribution: None,
          time_length: None,
        };
        (edge.key(), branch)
      })
      .collect(),
    backward: TimeBackward {
      subtree: graph.get_nodes().map(|node| (node.key(), None)).collect(),
      messages: graph.get_edges().map(|edge| (edge.key(), None)).collect(),
    },
    posterior: graph
      .get_nodes()
      .map(|node| (node.key(), NodePosterior::default()))
      .collect(),
  }
}

pub(crate) fn marginal_timetree_params() -> TimetreeParams {
  TimetreeParams {
    model: GtrModelName::JC69,
    alphabet_name: AlphabetName::Nuc,
    dense: None,
    gap_fill: GapFill::default(),
    branch_length_mode: BranchLengthMode::Marginal,
    no_indels: false,
    sequence_length: None,
    clock_rate: None,
    clock_std_dev: None,
    keep_root: true,
    reroot_spec: RerootSpec::default(),
    allow_negative_rate: false,
    clock_filter: 0.0,
    covariation: false,
    tip_slack: None,
    max_iter: 1,
    resolve_polytomies: false,
    keep_polytomies: false,
    relax: vec![],
    coalescent: None,
    coalescent_opt: false,
    coalescent_skyline: false,
    skyline_n_points: 0,
    skyline_stiffness: 0.0,
    coalescent_confidence: 0.0,
    gen_per_year: 0.0,
    n_branches_posterior: None,
    time_marginal: TimeMarginalMode::Never,
    confidence: false,
    include_leaves: false,
    impute_missing_data: false,
    report_ambiguous: false,
    zero_based: false,
    seed: None,
  }
}
