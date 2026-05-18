use crate::clock::clock_model::{ClockModel, ClockModelStats};
use crate::make_error;
use crate::partition::timetree::GraphTimetree;
use crate::partition::traits::PartitionTimetreeAll;
use crate::payload::timetree::EdgeTimetree;
use crate::payload::timetree::NodeTimetree;
use crate::payload::traits::TimetreeNode;
use crate::timetree::inference::runner::run_timetree;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::{info, warn};
use ordered_float::OrderedFloat;
use parking_lot::RwLock;
use serde::Serialize;
use statrs::function::erf::erf_inv;
use std::collections::BTreeMap;
use std::f64::consts::SQRT_2;
use std::path::Path;
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::{GraphNodeKey, Named, TimeConstraint};
use treetime_io::csv::CsvStructFileWriter;

/// Fraction of probability mass for the 90% confidence region.
/// Matches v0's `get_max_posterior_region(n, fraction=0.9)`.
const CI_FRACTION: f64 = 0.9;

/// Symmetric quantile bounds derived from CI_FRACTION, for equal-tailed CI.
/// Used by the rate contribution (which assumes Gaussian rate uncertainty)
/// and as fallback when HPD is unavailable.
/// v0 (clock_tree.py:1181): `((1-fraction)*0.5, 1.0-(1-fraction)*0.5)`
const CI_LOWER_QUANTILE: f64 = (1.0 - CI_FRACTION) * 0.5; // 0.05
const CI_UPPER_QUANTILE: f64 = 1.0 - (1.0 - CI_FRACTION) * 0.5; // 0.95

/// Quantify how clock rate uncertainty propagates to node date uncertainty.
///
/// Rate susceptibility analysis (Sagulenko, Puller & Neher 2018, Section 2.2):
/// node date uncertainty has two independent sources:
///
/// 1. **Mutation stochasticity** - captured by the marginal posterior from belief
///    propagation (backward/forward passes). Implemented in `extract_confidence_intervals`.
/// 2. **Clock rate uncertainty** - the regression slope has a standard error from
///    the Hessian inverse. Since times scale as `t = divergence / rate`, rate
///    uncertainty propagates to all node dates. Deeper nodes have higher sensitivity:
///    a 10% rate error shifts the root date by 10% of the tree depth.
///
/// This function addresses source (2) by re-running timetree inference three times:
/// at rate+sigma, rate-sigma, and central rate. Per-edge gamma values are scaled
/// proportionally (matching v0's `calc_rate_susceptibility` in `clock_tree.py:1010-1066`).
/// The resulting per-node date triples are stored sorted by date for later conversion
/// to confidence intervals by `date_uncertainty_due_to_rate`.
///
/// The third run (central rate) restores the graph to its pre-call state, so the
/// caller can proceed with further passes (e.g., final marginal reconstruction).
pub fn compute_rate_susceptibility(
  graph: &mut GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  clock_model: &ClockModel,
  coalescent_tc: Option<&Distribution>,
  rate_std: f64,
  no_indels: bool,
) -> Result<(), Report> {
  let current_rate = clock_model.clock_rate();

  // v0: upper = rate + rate_std, lower = max(0.1 * rate, rate - rate_std)
  // The 0.1*rate floor prevents sign reversal when rate_std > rate.
  let upper_rate = current_rate + rate_std;
  let lower_rate = (0.1 * current_rate).max(current_rate - rate_std);

  // Scale gammas by rate ratio, matching v0's gamma scaling approach.
  // v0 (clock_tree.py:1036-1039):
  //   n.branch_length_interpolator.gamma = n._orig_gamma * upper_rate / current_rate
  //
  // In run_timetree, each branch distribution is computed with:
  //   effective_clock_rate = clock_rate * gamma
  // So scaling gamma by (new_rate / current_rate) makes:
  //   effective_clock_rate = clock_rate * orig_gamma * new_rate / current_rate
  //                        = orig_gamma * new_rate
  // which is the desired effective rate at each branch.
  let original_gammas = save_gammas(graph);

  // Run 1: upper rate bound
  scale_gammas(graph, &original_gammas, upper_rate / current_rate);
  info!("Rate susceptibility: running with upper rate {upper_rate:.6e}");
  run_timetree(graph, partitions, clock_model, coalescent_tc, no_indels)
    .wrap_err("Rate susceptibility: timetree at upper rate failed")?;
  let upper_dates = collect_node_times(graph);

  // Run 2: lower rate bound
  scale_gammas(graph, &original_gammas, lower_rate / current_rate);
  info!("Rate susceptibility: running with lower rate {lower_rate:.6e}");
  run_timetree(graph, partitions, clock_model, coalescent_tc, no_indels)
    .wrap_err("Rate susceptibility: timetree at lower rate failed")?;
  let lower_dates = collect_node_times(graph);

  // Run 3: central rate (restores graph to pre-call state)
  scale_gammas(graph, &original_gammas, 1.0);
  info!("Rate susceptibility: running with central rate {current_rate:.6e}");
  run_timetree(graph, partitions, clock_model, coalescent_tc, no_indels)
    .wrap_err("Rate susceptibility: timetree at central rate failed")?;

  // Store sorted date triples per node.
  // v0 (clock_tree.py:1064): n.numdate_rate_variation.sort(key=lambda x: x[1])
  // The sort ensures [lower_date, central_date, upper_date] regardless of which
  // rate produced which date (deeper nodes may have inverted rate-date relationship).
  for node_ref in graph.get_nodes() {
    let node = node_ref.read_arc();
    let key = node.key();
    let mut payload = node.payload().write_arc();

    let central_date = payload.time();
    let upper_date = upper_dates.get(&key).copied();
    let lower_date = lower_dates.get(&key).copied();

    if let (Some(c), Some(u), Some(l)) = (central_date, upper_date, lower_date) {
      let mut dates = [l, c, u];
      dates.sort_by_key(|x| OrderedFloat(*x));
      payload.rate_susceptibility_dates = Some(dates);
    }
  }

  info!("Rate susceptibility analysis completed");
  Ok(())
}

/// Convert per-node rate-variation date triples to a confidence interval.
///
/// The date triple `[lower, central, upper]` represents node dates inferred at
/// clock_rate - rate_std, clock_rate, and clock_rate + rate_std respectively
/// (after sorting by date value).
///
/// The mapping from +/-1 sigma rate variation to arbitrary quantile bounds uses
/// the probit function (inverse normal CDF), matching v0's
/// `date_uncertainty_due_to_rate` (clock_tree.py:1068-1088):
///
/// ```text
/// z = sqrt(2) * erf_inv(2p - 1)     // probit: quantile p to z-score
/// ci_lower = central + z(p_lo) * |lower - central|
/// ci_upper = central + z(p_hi) * |upper - central|
/// ```
///
/// With CI_LOWER_QUANTILE=0.05 and CI_UPPER_QUANTILE=0.95 (90% CI):
/// z(0.05)=-1.645, z(0.95)=+1.645. The +-1 sigma date variation (from
/// rate +/- 1 std) is scaled to the 90% CI level: the date deviation at
/// +/-1 sigma is multiplied by 1.645 to get the 5/95 percentile bounds.
/// This assumes the rate-to-date mapping is approximately linear (valid
/// when rate_std << rate) and that the rate estimate is approximately
/// Gaussian (from the regression CLT).
pub(crate) fn date_uncertainty_due_to_rate(dates: [f64; 3], interval: (f64, f64)) -> (f64, f64) {
  let [lower, central, upper] = dates;
  let z_lower = quantile_to_zscore(interval.0);
  let z_upper = quantile_to_zscore(interval.1);
  // v0 (clock_tree.py:1085):
  //   return np.array([c + x * np.abs(y - c) for x, y in zip(nsig, (l, u))])
  let ci_lower = central + z_lower * (lower - central).abs();
  let ci_upper = central + z_upper * (upper - central).abs();
  (ci_lower, ci_upper)
}

/// Extract confidence intervals combining marginal posterior and rate susceptibility.
///
/// When both marginal distributions and rate susceptibility data are present,
/// combines the two independent uncertainty sources via quadrature sum
/// (Sagulenko et al. 2018, Section 2.2). This treats mutation stochasticity
/// and clock rate uncertainty as independent Gaussian-like contributions:
///
/// ```text
/// combined_lower = center - sqrt((rate_lo - center)^2 + (marginal_lo - center)^2)
/// combined_upper = center + sqrt((rate_hi - center)^2 + (marginal_hi - center)^2)
/// ```
///
/// When only one source is present, uses that source alone.
/// When neither is present, returns point interval [date, date].
pub fn extract_confidence_intervals(graph: &GraphTimetree) -> Vec<NodeConfidenceInterval> {
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node_ref| {
      let node = node_ref.read_arc();
      let key = node.key();
      let payload = node.payload().read_arc();
      let name = payload.name().map_or_else(String::new, |n| n.as_ref().to_owned());
      let date = payload.time?;

      // Source 1: mutation stochasticity from marginal posterior HPD region.
      // v0 uses get_max_posterior_region(fraction=0.9): highest posterior density
      // region, the narrowest interval containing 90% of the probability mass.
      // For symmetric distributions this equals equal-tailed CI; for skewed
      // distributions (nodes near tree boundaries) HPD is narrower and
      // centered on the peak.
      let mutation_contribution = payload
        .time_distribution()
        .as_ref()
        .and_then(|dist| dist.hpd_region(CI_FRACTION));

      // Source 2: clock rate uncertainty from rate susceptibility analysis
      let rate_contribution = payload
        .rate_susceptibility_dates
        .map(|dates| date_uncertainty_due_to_rate(dates, (CI_LOWER_QUANTILE, CI_UPPER_QUANTILE)));

      let (lower, upper) = if rate_contribution.is_none() && mutation_contribution.is_none() {
        // No CI data available: return point estimate
        (date, date)
      } else {
        // Physical limits from distribution domain (or unbounded if no distribution)
        let limits = payload
          .time_distribution()
          .as_ref()
          .map_or((f64::NEG_INFINITY, f64::INFINITY), |dist| dist.time_bounds());
        combine_confidence(date, limits, rate_contribution, mutation_contribution)
      };

      // Postcondition: point estimate must lie within confidence interval.
      // When only one contribution is present, combine_confidence returns its
      // raw bounds, which are centered on that contribution's own reference
      // point (rate susceptibility central date or HPD peak), not on `date`.
      // If those differ after the final marginal pass, the raw CI may not
      // bracket the point estimate. Clamp to ensure validity.
      let lower = lower.min(date);
      let upper = upper.max(date);

      Some(NodeConfidenceInterval {
        key,
        name,
        date,
        lower,
        upper,
      })
    })
    .sorted_by_key(|ci| ci.key)
    .collect_vec()
}

/// Confidence interval for a node's inferred date.
#[derive(Debug, Clone, Serialize)]
pub struct NodeConfidenceInterval {
  /// Graph node key for stable lookup across naming schemes.
  #[serde(skip)]
  pub key: GraphNodeKey,
  pub name: String,
  pub date: f64,
  pub lower: f64,
  pub upper: f64,
}

/// Write confidence intervals for node dates to TSV file.
pub fn write_confidence_intervals(intervals: &[NodeConfidenceInterval], filepath: &Path) -> Result<(), Report> {
  let mut writer = CsvStructFileWriter::new(filepath, b'\t')?;
  intervals.iter().try_for_each(|ci| writer.write(ci))
}

/// Combine two confidence interval contributions via quadrature sum.
///
/// Given a center date and physical limits, combines up to two independent
/// confidence interval contributions (e.g., from mutation uncertainty and rate uncertainty).
/// Returns bounds clipped to the physical limits.
///
/// v0: `combine_confidence()` in `clock_tree.py:1090-1101`.
///
/// Formula:
/// - `lower = center - sqrt((c1_lo - center)^2 + (c2_lo - center)^2)`
/// - `upper = center + sqrt((c1_hi - center)^2 + (c2_hi - center)^2)`
pub(crate) fn combine_confidence(
  center: f64,
  limits: (f64, f64),
  c1: Option<(f64, f64)>,
  c2: Option<(f64, f64)>,
) -> (f64, f64) {
  let (min_val, max_val) = match (c1, c2) {
    (None, None) => return limits,
    (Some(c), None) | (None, Some(c)) => c,
    (Some(c1), Some(c2)) => {
      let min_val = center - (c1.0 - center).hypot(c2.0 - center);
      let max_val = center + (c1.1 - center).hypot(c2.1 - center);
      (min_val, max_val)
    },
  };

  (limits.0.max(min_val), limits.1.min(max_val))
}

/// Determine rate standard deviation for rate susceptibility analysis.
///
/// Priority: `--clock-std-dev` (user-specified) > regression covariance matrix.
///
/// v0 gating (wrappers.py:453-469, clock_tree.py:1020-1026):
/// - `--clock-std-dev` provides the value directly
/// - `--covariation` enables `valid_confidence=True`, allowing `sqrt(cov[0,0])`
/// - Without either, v0 disables confidence estimation entirely
///
/// Returns `None` when rate susceptibility should not run.
pub(crate) fn determine_rate_std(
  clock_std_dev: Option<f64>,
  covariation: bool,
  clock_model: &ClockModel,
) -> Result<Option<f64>, Report> {
  if let Some(std_dev) = clock_std_dev {
    if std_dev <= 0.0 {
      return make_error!("--clock-std-dev must be positive, got {std_dev}");
    }
    return Ok(Some(std_dev));
  }

  if !covariation {
    return Ok(None);
  }

  // Derive rate_std from regression covariance matrix.
  // v0 (clock_tree.py:1026): rate_std = np.sqrt(self.clock_model['cov'][0, 0])
  // The cov[0,0] entry is the variance of the slope (clock rate) from GLS regression.
  // Only available with --covariation (phylogenetic GLS), which sets valid_confidence=True.
  match clock_model.stats() {
    ClockModelStats::Estimated(stats) => {
      let rate_variance = stats.cov[[0, 0]];
      if rate_variance <= 0.0 {
        warn!(
          "Rate variance from regression covariance is non-positive ({rate_variance:.4e}), skipping rate susceptibility"
        );
        return Ok(None);
      }
      Ok(Some(rate_variance.sqrt()))
    },
    ClockModelStats::Fixed => Ok(None),
  }
}

/// Convert a probability quantile to a z-score via the probit function.
///
/// `z = sqrt(2) * erf_inv(2p - 1)` is the quantile function (probit) of the
/// standard normal distribution. Maps p in (0, 1) to the z-score such that
/// Phi(z) = p, where Phi is the standard normal CDF.
///
/// Returns 0.0 for boundary values p=0 or p=1, matching v0's guard:
/// `if x * (1.0 - x) else 0` (clock_tree.py:1083).
pub(crate) fn quantile_to_zscore(p: f64) -> f64 {
  // v0 guard: if x * (1.0 - x) is falsy (i.e. x=0 or x=1), return 0
  if p * (1.0 - p) == 0.0 {
    return 0.0;
  }
  SQRT_2 * erf_inv(2.0 * p - 1.0)
}

/// Save per-edge gamma values for later restoration.
fn save_gammas(graph: &GraphTimetree) -> Vec<(GraphEdgeKey, f64)> {
  graph
    .get_edges()
    .into_iter()
    .map(|edge_ref| {
      let edge = edge_ref.read_arc();
      (edge.key(), edge.payload().read_arc().gamma)
    })
    .collect_vec()
}

/// Scale all edge gammas: `new_gamma = original_gamma * scale_factor`.
fn scale_gammas(graph: &GraphTimetree, original_gammas: &[(GraphEdgeKey, f64)], scale_factor: f64) {
  for &(key, orig_gamma) in original_gammas {
    if let Some(edge_ref) = graph.get_edge(key) {
      edge_ref.write_arc().payload().write_arc().gamma = orig_gamma * scale_factor;
    }
  }
}

/// Collect per-node time estimates keyed by node key.
fn collect_node_times(graph: &GraphTimetree) -> BTreeMap<GraphNodeKey, f64> {
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node_ref| {
      let node = node_ref.read_arc();
      Some((node.key(), node.payload().read_arc().time()?))
    })
    .collect()
}
