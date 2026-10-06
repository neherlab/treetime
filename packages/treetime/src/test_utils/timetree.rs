use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::RerootSpec;
use crate::coalescent::node_time::{CoalescentNodeTime, CoalescentNodeTimes};
use crate::gtr::get_gtr::GtrModelName;
use crate::optimize::params::BranchLengthMode;
use crate::progress::{LogEvent, LogLevel, LogSink};
use crate::test_utils::find_node_key_by_name;
use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
use crate::timetree::inference::result::{BranchLikelihood, NodePosterior, TimeInference, given_times};
use crate::timetree::params::TimeMarginalMode;
use crate::timetree::params::TimetreeParams;
use eyre::Report;
use parking_lot::Mutex;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::MaxGridPoints;

pub(crate) fn constraint_coalescent_node_times(
  graph: &Graph,
  constraints: &DateConstraints,
) -> Result<CoalescentNodeTimes, Report> {
  let bad_branches = derive_bad_branches(graph, constraints, &bad_leaves(graph, constraints, &BTreeSet::new()))?;
  let times = given_times(graph, constraints)?
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
    branches: unknown_branches(graph),
    posterior: graph
      .get_nodes()
      .map(|node| (node.key(), NodePosterior::default()))
      .collect(),
  }
}

pub(crate) fn marginal_timetree_params() -> TimetreeParams {
  TimetreeParams {
    model: GtrModelName::JC69,
    dense: None,
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
    report_ambiguous: true,
    impute_missing_data: false,
    sequence_outputs_requested: false,
    seed: 0,
    max_grid_points: MaxGridPoints::default(),
  }
}

pub(crate) fn unknown_branches(graph: &Graph) -> BTreeMap<GraphEdgeKey, BranchLikelihood> {
  graph
    .get_edges()
    .map(|edge| {
      let branch = BranchLikelihood {
        distribution: None,
        time_length: None,
      };
      (edge.key(), branch)
    })
    .collect()
}

pub(crate) fn parent_edge_key(graph: &Graph, target_key: GraphNodeKey) -> GraphEdgeKey {
  graph
    .get_edges()
    .find(|edge| edge.target() == target_key)
    .expect("node must have a parent edge")
    .key()
}

pub(crate) fn point_date_constraints(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  dates: &[(&str, f64)],
) -> DateConstraints {
  let by_node = dates
    .iter()
    .map(|(name, date)| {
      let key = find_node_key_by_name(graph, names, name).expect("dated node must exist");
      (key, Some(Arc::new(Distribution::point(*date, 0.0))))
    })
    .collect();
  DateConstraints { by_node }
}

#[derive(Default)]
pub(crate) struct RecordingLog {
  events: Mutex<Vec<LogEvent>>,
}

impl RecordingLog {
  pub(crate) fn warnings(&self) -> Vec<String> {
    self
      .events
      .lock()
      .iter()
      .filter(|event| event.level == LogLevel::Warn)
      .map(|event| event.message.clone())
      .collect()
  }
}

impl LogSink for RecordingLog {
  fn log(&self, level: LogLevel, message: &str) {
    self.events.lock().push(LogEvent {
      level,
      message: message.to_owned(),
    });
  }

  fn log_enabled(&self, _level: LogLevel) -> bool {
    true
  }
}
