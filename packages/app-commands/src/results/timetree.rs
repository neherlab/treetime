use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
use crate::json_float::JsonFloat;
use crate::results::clock::{RootToTip, is_fixed, root_to_tip};
use crate::results::tree::{DateInterval, ResultTree};
use crate::results::year_date::YearDate;
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::rtt::ClockRegressionResult;
use treetime::timetree::coalescent::CoalescentSegmentRow;
use treetime::timetree::coalescent_timescale::{CoalescentMode, coalescent_mode};
use treetime::timetree::convergence::metrics::ConvergenceMetrics;
use treetime_schema::skip_serializing_optionals;
use util_augur_node_data_json::AugurNodeDataJsonClock;

const INTERVAL_EDGE_FRACTION: f64 = 0.05;

/// Results of a `timetree` run.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct TimetreeResults {
  /// Estimates of the time tree; absent when the run wrote no Auspice tree.
  pub estimates: Option<TimetreeEstimates>,
  /// Samples and line of the final clock model; absent when the run wrote no clock regression table.
  pub root_to_tip: Option<RootToTip>,
  /// Convergence values of every iteration, from the tracelog.
  pub iterations: Vec<IterationRow>,
  /// Segments of the coalescent time scale, from the coalescent table.
  pub skyline: Vec<SkylineSegment>,
}

/// Estimates of a `timetree` run.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct TimetreeEstimates {
  /// Date of the root.
  pub root_date: Option<YearDate>,
  /// Confidence interval of the root date.
  pub root_interval: Option<DateInterval>,
  /// Whether the root date lies within 5% of the interval width from a bound of its interval.
  pub root_near_interval_edge: bool,
  /// Clock rate in substitutions per site per year.
  pub clock_rate: Option<f64>,
  /// Standard deviation of the clock rate: the estimate's, or the one given with a fixed rate.
  pub clock_rate_std: Option<f64>,
  /// Whether the clock rate was fixed by the user instead of estimated.
  pub clock_rate_fixed: bool,
  /// Correlation coefficient of the final clock model; absent for a fixed rate.
  pub r: Option<f64>,
  /// Coefficient of determination of the final clock model.
  pub r_squared: Option<f64>,
  /// Number of samples in the tree.
  pub samples: usize,
  /// Number of samples the clock model left out: without a usable date or clock outliers.
  pub excluded_samples: usize,
  /// Coalescent prior the run used.
  pub coalescent_prior: CoalescentPrior,
  /// Relaxed clock the run used.
  pub relaxed_clock: Option<RelaxedClock>,
  /// Total log likelihood of the last iteration.
  pub log_likelihood: Option<JsonFloat>,
  /// Number of iterations the tracelog records.
  pub iterations: usize,
}

/// Coalescent prior of a `timetree` run.
#[derive(Clone, Copy, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
#[schemars(tag = "kind", rename_all = "kebab-case")]
#[deser(tag = "kind", rename_all = "kebab-case")]
pub enum CoalescentPrior {
  /// No coalescent prior.
  None,
  /// Constant population size with a fixed time scale.
  Fixed {
    /// Coalescent time scale in years.
    tc: f64,
  },
  /// Constant population size with an optimized time scale.
  Optimized,
  /// Piecewise-constant population size.
  Skyline {
    /// Number of grid points.
    points: usize,
    /// Stiffness of the skyline.
    stiffness: f64,
  },
}

/// Parameters of a relaxed clock.
#[derive(Clone, Copy, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
pub struct RelaxedClock {
  /// Slack: how far the rate of a branch may vary.
  pub slack: f64,
  /// Coupling: how strongly the rates of parent and child branches are tied.
  pub coupling: f64,
}

/// Convergence values of one iteration.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct IterationRow {
  /// Iteration number, from 0.
  pub iteration: usize,
  /// Largest change of a node time, in years.
  pub max_time_change: Option<JsonFloat>,
  /// Root-mean-square change of the node times, in years.
  pub rms_time_change: Option<JsonFloat>,
  /// Log likelihood of the sequences.
  pub log_lh_seq: Option<JsonFloat>,
  /// Log likelihood of the node positions.
  pub log_lh_pos: Option<JsonFloat>,
  /// Log likelihood of the coalescent prior.
  pub log_lh_coal: Option<JsonFloat>,
  /// Total log likelihood.
  pub log_lh_total: Option<JsonFloat>,
}

/// One segment of the coalescent time scale.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
pub struct SkylineSegment {
  /// Start of the segment, as a decimal year.
  pub start: f64,
  /// End of the segment, as a decimal year.
  pub end: f64,
  /// Coalescent time scale in years.
  pub tc: Band,
  /// Effective population size.
  pub ne: Band,
}

/// An estimate with an optional confidence band.
#[derive(Clone, Copy, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct Band {
  /// Point estimate.
  pub value: f64,
  /// Lower bound of the band.
  pub lower: Option<f64>,
  /// Upper bound of the band.
  pub upper: Option<f64>,
}

pub struct TimetreeOutputs<'a> {
  pub tree: Option<&'a ResultTree>,
  pub clock_model: Option<&'a ClockModel>,
  pub clock_rows: Option<&'a [ClockRegressionResult]>,
  pub node_data_clock: Option<&'a AugurNodeDataJsonClock>,
  pub trace: &'a [ConvergenceMetrics],
  pub coalescent: &'a [CoalescentSegmentRow],
}

pub fn timetree_results(outputs: &TimetreeOutputs<'_>, config: &TreetimeTimetreeArgsRaw) -> TimetreeResults {
  TimetreeResults {
    estimates: outputs.tree.map(|tree| timetree_estimates(tree, outputs, config)),
    root_to_tip: outputs
      .clock_rows
      .map(|rows| root_to_tip(outputs.tree, outputs.clock_model, rows)),
    iterations: outputs
      .trace
      .iter()
      .enumerate()
      .map(|(iteration, metrics)| IterationRow {
        iteration,
        max_time_change: metrics.max_time_change.map(JsonFloat),
        rms_time_change: metrics.rms_time_change.map(JsonFloat),
        log_lh_seq: metrics.log_lh_seq.map(JsonFloat::from),
        log_lh_pos: metrics.log_lh_pos.map(JsonFloat::from),
        log_lh_coal: metrics.log_lh_coal.map(JsonFloat::from),
        log_lh_total: metrics.log_lh_total.map(JsonFloat::from),
      })
      .collect(),
    skyline: outputs
      .coalescent
      .iter()
      .map(|row| SkylineSegment {
        start: row.segment_start,
        end: row.segment_end,
        tc: Band {
          value: row.tc_value,
          lower: row.tc_lower,
          upper: row.tc_upper,
        },
        ne: Band {
          value: row.ne_value,
          lower: row.ne_lower,
          upper: row.ne_upper,
        },
      })
      .collect(),
  }
}

pub fn coalescent_prior(config: &TreetimeTimetreeArgsRaw) -> CoalescentPrior {
  match coalescent_mode(config.coalescent, config.coalescent_opt, config.coalescent_skyline) {
    CoalescentMode::Skyline => CoalescentPrior::Skyline {
      points: config.skyline_n_points,
      stiffness: config.skyline_stiffness,
    },
    CoalescentMode::Constant => CoalescentPrior::Optimized,
    CoalescentMode::Fixed(tc) => CoalescentPrior::Fixed { tc },
    CoalescentMode::Disabled => CoalescentPrior::None,
  }
}

pub fn relaxed_clock(config: &TreetimeTimetreeArgsRaw) -> Option<RelaxedClock> {
  match config.relax.as_slice() {
    &[slack, coupling] => Some(RelaxedClock { slack, coupling }),
    _ => None,
  }
}

fn timetree_estimates(
  tree: &ResultTree,
  outputs: &TimetreeOutputs<'_>,
  config: &TreetimeTimetreeArgsRaw,
) -> TimetreeEstimates {
  let root = tree.root();
  let fixed = outputs.clock_model.is_some_and(is_fixed);
  let r = outputs.clock_model.and_then(ClockModel::r_val);
  TimetreeEstimates {
    root_date: root.date.clone(),
    root_interval: root.date_interval.clone(),
    root_near_interval_edge: root
      .date
      .as_ref()
      .zip(root.date_interval.as_ref())
      .is_some_and(|(date, interval)| near_interval_edge(date.year, interval)),
    clock_rate: outputs
      .clock_model
      .map(ClockModel::clock_rate)
      .or_else(|| outputs.node_data_clock.map(|clock| clock.rate)),
    clock_rate_std: if fixed {
      config.clock_std_dev
    } else {
      outputs.node_data_clock.and_then(|clock| clock.rate_std)
    },
    clock_rate_fixed: fixed,
    r,
    r_squared: outputs.clock_model.and_then(ClockModel::r_squared),
    samples: tree.tips().count(),
    excluded_samples: tree.tips().filter(|tip| tip.excluded == Some(true)).count(),
    coalescent_prior: coalescent_prior(config),
    relaxed_clock: relaxed_clock(config),
    log_likelihood: outputs
      .trace
      .last()
      .and_then(|metrics| metrics.log_lh_total.map(JsonFloat::from)),
    iterations: outputs.trace.len(),
  }
}

fn near_interval_edge(date: f64, interval: &DateInterval) -> bool {
  let (lower, upper) = (interval.lower.year, interval.upper.year);
  (date - lower).min(upper - date) <= INTERVAL_EDGE_FRACTION * (upper - lower)
}
