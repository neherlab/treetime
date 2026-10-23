use crate::clock::clock_model::ClockModel;
use crate::clock::date_constraints::DateConstraints;
use crate::coalescent::mean_field::CoalescentPrior;
use crate::optimize::indel::IndelPolicy;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::reconstruction::DenseReconstruction;
use crate::partition::set::PartitionSet;
use crate::progress::NoopProgress;
use crate::test_utils::{
  BranchTime, CaseSpec, NodeFactors, ReferenceInput, ReferenceLogZ, SimulatedCase, gtr_from_model,
  point_coalescent_prior_by_key, reference_log_z, simulated_sequence_source,
};
use crate::timetree::branch_model::BranchModel;
use crate::timetree::inference::bad_branches::bad_leaves;
use crate::timetree::inference::result::{BranchLikelihood, PosteriorSummaries};
use crate::timetree::inference::runner::{TimeInferenceInputs, run_timetree};
use crate::timetree::optimization::marginal_clock_rate::{
  MarginalRateInputs, RateWindow, branches_at_rate, marginal_log_lh_at_rate,
};
use crate::timetree::optimization::relaxed_clock::unit_gammas;
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_grid::MaxGridPoints;
use treetime_simulation::branch_lengths::realized;
use treetime_simulation::tc_schedule::TcSchedule;
use treetime_utils::make_error;

pub(crate) struct FrozenRateState {
  pub spec: CaseSpec,
  pub case: SimulatedCase,
  pub tc: f64,
  pub tip_dates: BTreeMap<GraphNodeKey, f64>,
  pub node_dates: BTreeMap<GraphNodeKey, f64>,
  pub constraints: DateConstraints,
  pub bad_branches: BTreeMap<GraphNodeKey, bool>,
  pub branches: BTreeMap<GraphEdgeKey, BranchLikelihood>,
  pub prior: CoalescentPrior,
}

impl FrozenRateState {
  pub(crate) fn new(spec: &CaseSpec) -> Result<Self, Report> {
    let tc = spec.schedule.tc(0.0);
    if spec.schedule != TcSchedule::constant(tc)? {
      return make_error!(
        "The frozen rate state needs a constant Tc schedule, got {:?}",
        spec.schedule
      );
    }
    let case = spec.simulate()?;
    let graph = &case.tree.graph;
    let node_dates: BTreeMap<GraphNodeKey, f64> = graph
      .get_nodes()
      .map(|node| (node.key(), case.tree.date(node.key())))
      .collect();
    let tip_dates: BTreeMap<GraphNodeKey, f64> = graph
      .get_leaves()
      .map(|leaf| (leaf.key(), node_dates[&leaf.key()]))
      .collect();
    let constraints = DateConstraints {
      by_node: tip_dates
        .iter()
        .map(|(&key, &date)| (key, Some(Arc::new(Distribution::point(date, 0.0)))))
        .collect(),
    };
    let (bad_branches, branches) = true_rate_branches(spec, &case, &constraints)?;
    let prior = point_coalescent_prior_by_key(graph, &node_dates, tc)?;
    Ok(Self {
      spec: spec.clone(),
      case,
      tc,
      tip_dates,
      node_dates,
      constraints,
      bad_branches,
      branches,
      prior,
    })
  }

  pub(crate) fn profile(&self, rates: &[f64], with_prior: bool, spacing: f64) -> Result<Vec<RateProfilePoint>, Report> {
    rates
      .iter()
      .map(|&rate| {
        self
          .point(rate, with_prior, spacing)
          .wrap_err_with(|| format!("When profiling log Z at clock rate {rate:e}"))
      })
      .collect()
  }

  /// TreeTime's log Z at `rate` under the frozen prior, without the reference
  pub(crate) fn log_z(&self, rate: f64, with_prior: bool) -> Result<Option<f64>, Report> {
    let inputs = MarginalRateInputs {
      graph: &self.case.tree.graph,
      date_constraints: &self.constraints,
      coalescent: with_prior.then_some(&self.prior),
      bad_branches: &self.bad_branches,
      branches: &self.branches,
      reference_rate: self.spec.mu,
      window: RateWindow::InteriorMaximum,
      max_points: MaxGridPoints::default(),
    };
    marginal_log_lh_at_rate(&inputs, rate)
  }

  pub(crate) fn point(&self, rate: f64, with_prior: bool, spacing: f64) -> Result<RateProfilePoint, Report> {
    let prior = with_prior.then_some(&self.prior);
    let treetime_log_z = self.log_z(rate, with_prior)?;
    let log_normalization = prior.map_or(0.0, prior_log_normalization);
    let branches = BranchTime::from_likelihoods(&branches_at_rate(&self.branches, self.spec.mu, rate)?)?;
    let reference = reference_log_z(&ReferenceInput {
      graph: &self.case.tree.graph,
      tip_dates: &self.tip_dates,
      branches: &branches,
      factors: prior.map_or(NodeFactors::None, NodeFactors::Prior),
      spacing,
    })?;
    Ok(RateProfilePoint {
      rate,
      treetime_log_z,
      treetime_raw_log_z: treetime_log_z.map(|log_z| log_z - log_normalization),
      log_normalization,
      reference,
    })
  }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct RateProfilePoint {
  pub rate: f64,
  pub treetime_log_z: Option<f64>,
  pub treetime_raw_log_z: Option<f64>,
  pub log_normalization: f64,
  pub reference: ReferenceLogZ,
}

fn prior_log_normalization(prior: &CoalescentPrior) -> f64 {
  prior.log_normalization()
}

fn true_rate_branches(
  spec: &CaseSpec,
  case: &SimulatedCase,
  constraints: &DateConstraints,
) -> Result<(BTreeMap<GraphNodeKey, bool>, BTreeMap<GraphEdgeKey, BranchLikelihood>), Report> {
  let graph = &case.tree.graph;
  let source = simulated_sequence_source(graph, &case.tree.names, &case.sequences.sequences)?;
  let partition = PartitionMarginalDense::new(source.alphabet, graph, &source.nodes)?;
  let gtr = gtr_from_model(&case.model, 1.0)?;
  let lengths = realized(graph, &case.sequences.sequences)?;
  let (partitions, _) =
    PartitionSet::Dense(vec![DenseReconstruction::seeded(partition, gtr)]).marginal_update(graph, &lengths)?;
  let branch_lengths = lengths.iter().map(|(&edge, &length)| (edge, Some(length))).collect();
  let inference = run_timetree(
    &TimeInferenceInputs {
      graph,
      date_constraints: constraints,
      leaf_bad_branches: &bad_leaves(graph, constraints, &BTreeSet::new()),
      gammas: &unit_gammas(graph),
      branch_model: &BranchModel::Marginal(partitions),
      branch_lengths: &branch_lengths,
      names: &case.tree.names,
      clock_model: &ClockModel::for_testing(spec.mu, 0.0),
      clock_rate_fixed: true,
      indel_policy: IndelPolicy::Include,
      max_grid_points: MaxGridPoints::default(),
      posterior_summaries: PosteriorSummaries::default(),
    },
    None,
    &NoopProgress,
  )?;
  Ok((inference.bad_branches, inference.branches))
}
