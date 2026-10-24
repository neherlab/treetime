use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::RerootSpec;
use crate::coalescent::mean_field::{CoalescentPrior, LineageExpectation};
use crate::gtr::get_gtr::{GtrModel, JC69Params};
use crate::optimize::indel::IndelPolicy;
use crate::optimize::partition_rates::PartitionRates;
use crate::partition::split::SequenceInput;
use crate::seq::alignment::SequenceSource;
use crate::test_utils::find_node_key_by_name;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::coalescent_timescale::CoalescentMode;
use crate::timetree::inference::input_branch_length::InputBranchModel;
use crate::timetree::inference::result::{
  BranchLikelihood, NodePosterior, PosteriorSummaries, PosteriorView, TimeDistribution, TimeInference, given_times,
};
use crate::timetree::params::{BranchLengthMode, TimetreeParams};
use crate::timetree::pipeline::TimetreeSequences;
use eyre::Report;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_grid::MaxGridPoints;

pub(crate) const JC69_NUCLEOTIDE_SATURATED_DIFFERENCE: f64 = 0.75;

pub(crate) fn jc69_input_branch_model(sequence_length: usize) -> BranchModel {
  BranchModel::Input(InputBranchModel::new(
    sequence_length,
    JC69_NUCLEOTIDE_SATURATED_DIFFERENCE,
  ))
}

pub(crate) fn jc69_timetree_sequences(sources: Vec<SequenceSource>) -> Result<TimetreeSequences, Report> {
  Ok(TimetreeSequences::Alignments(SequenceInput::whole_alignments(
    sources,
    &GtrModel::JC69(JC69Params::default()),
  )?))
}

pub(crate) fn constraint_expectation(
  graph: &Graph,
  constraints: &DateConstraints,
) -> Result<Arc<LineageExpectation>, Report> {
  let posterior = given_times(graph, constraints)?
    .into_iter()
    .map(|(key, time)| {
      let entry = NodePosterior {
        likely_time: time,
        time,
        ..NodePosterior::default()
      };
      (key, entry)
    })
    .collect();
  LineageExpectation::from_posterior(graph, &posterior)
}

pub(crate) fn empty_time_inference(graph: &Graph) -> TimeInference {
  TimeInference {
    log_lh_time: None,
    bad_branches: graph.get_nodes().map(|node| (node.key(), false)).collect(),
    branches: unknown_branches(graph),
    posterior: graph
      .get_nodes()
      .map(|node| (node.key(), NodePosterior::default()))
      .collect(),
    posterior_summaries: PosteriorSummaries::default(),
  }
}

pub(crate) fn marginal_timetree_params() -> TimetreeParams {
  TimetreeParams {
    dense: None,
    branch_length_mode: BranchLengthMode::Marginal,
    indel_policy: IndelPolicy::Include,
    partition_rates: PartitionRates::Estimate,
    clock_rate: None,
    clock_std_dev: None,
    keep_root: true,
    reroot_spec: RerootSpec::default(),
    allow_negative_rate: false,
    clock_filter: 0.0,
    max_iter: 1,
    resolve_polytomies: false,
    relax: vec![],
    coalescent: CoalescentMode::Disabled,
    n_skyline: 0,
    skyline_stiffness: 0.0,
    coalescent_confidence: 0.0,
    gen_per_year: 0.0,
    confidence: false,
    include_leaves: false,
    report_ambiguous: true,
    impute_missing_data: false,
    sequence_outputs_requested: false,
    seed: 0,
    max_grid_points: MaxGridPoints::default(),
  }
}

pub(crate) fn full_posterior(_: GraphNodeKey, view: PosteriorView<'_>) -> Result<FullPosterior, Report> {
  Ok(FullPosterior {
    distribution: view.distribution.cloned(),
    likely_time: view.likely_time,
    time: view.time,
    contradicted: view.contradicted,
    subtree_support: view.subtree_support,
  })
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FullPosterior {
  pub distribution: TimeDistribution,
  pub likely_time: Option<f64>,
  pub time: Option<f64>,
  pub contradicted: bool,
  pub subtree_support: Option<(f64, f64)>,
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

pub(crate) fn point_date_constraints(graph: &Graph, names: &NodeNames, dates: &[(&str, f64)]) -> DateConstraints {
  let by_node = dates
    .iter()
    .map(|(name, date)| {
      let key = find_node_key_by_name(graph, names, name).expect("dated node must exist");
      (key, Some(Arc::new(Distribution::point(*date, 0.0))))
    })
    .collect();
  DateConstraints { by_node }
}

pub(crate) fn point_coalescent_prior(
  graph: &Graph,
  names: &NodeNames,
  times: &[(&str, f64)],
  tc: f64,
) -> Result<CoalescentPrior, Report> {
  let times = times
    .iter()
    .map(|(name, time)| {
      let key = find_node_key_by_name(graph, names, name).expect("timed node must exist");
      (key, *time)
    })
    .collect();
  point_coalescent_prior_by_key(graph, &times, tc)
}

pub(crate) fn point_coalescent_prior_by_key(
  graph: &Graph,
  times: &BTreeMap<GraphNodeKey, f64>,
  tc: f64,
) -> Result<CoalescentPrior, Report> {
  let posterior = times
    .iter()
    .map(|(&key, &time)| {
      let entry = NodePosterior {
        likely_time: Some(time),
        time: Some(time),
        ..NodePosterior::default()
      };
      (key, entry)
    })
    .collect();
  let expectation = LineageExpectation::from_posterior(graph, &posterior)?;
  CoalescentPrior::new(&expectation, &Distribution::constant(tc))
}
