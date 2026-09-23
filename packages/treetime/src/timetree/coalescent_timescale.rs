use crate::coalescent::lineage_counts::compute_lineage_counts;
use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
use crate::make_error;
use crate::timetree::coalescent::{
  CoalescentBand, CoalescentInputs, CoalescentOutput, CoalescentOutputMode, CoalescentSolve,
};
use eyre::{Report, WrapErr};
use ndarray::{Array1, array};
use treetime_distribution::Distribution;
use treetime_graph::graph::Graph;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::make_report;

pub fn coalescent_mode(coalescent: Option<f64>, coalescent_opt: bool, coalescent_skyline: bool) -> CoalescentMode {
  if coalescent_skyline {
    CoalescentMode::Skyline
  } else if coalescent_opt {
    CoalescentMode::Constant
  } else if let Some(tc) = coalescent {
    CoalescentMode::Fixed(tc)
  } else {
    CoalescentMode::Disabled
  }
}

pub fn coalescent_timescale(
  mode: CoalescentMode,
  graph: &Graph,
  skyline_params: &SkylineParams,
  node_times: &CoalescentNodeTimes,
) -> Result<CoalescentTimescale, Report> {
  let mode = match mode {
    CoalescentMode::Disabled => CoalescentMode::Constant,
    mode @ (CoalescentMode::Fixed(_) | CoalescentMode::Constant | CoalescentMode::Skyline) => mode,
  };
  estimate_coalescent_tc(mode, graph, skyline_params, node_times)
    .wrap_err("Failed to estimate the coalescent timescale")?
    .ok_or_else(|| make_report!("A coalescent Tc is required, but {mode:?} yielded none"))
}

pub fn estimate_coalescent_tc(
  mode: CoalescentMode,
  graph: &Graph,
  skyline_params: &SkylineParams,
  node_times: &CoalescentNodeTimes,
) -> Result<Option<CoalescentTimescale>, Report> {
  let n_points = match mode {
    CoalescentMode::Disabled => return Ok(None),
    CoalescentMode::Fixed(tc) => return fixed_timescale(tc, graph, node_times).map(Some),
    CoalescentMode::Constant => 1,
    CoalescentMode::Skyline => skyline_params.n_points,
  };
  let result = optimize_skyline(
    graph,
    &SkylineParams {
      n_points,
      ..skyline_params.clone()
    },
    node_times,
  )?;
  Ok(Some(CoalescentTimescale {
    distribution: result.tc_distribution,
    schedule: result.tc_schedule,
    report: Some(CoalescentTcReport {
      segment_boundaries: result.segment_boundaries,
      band: Some(CoalescentReportBand {
        lower: result.tc_lower_bounds,
        upper: result.tc_upper_bounds,
      }),
      log_likelihood: Some(result.log_likelihood.value()),
    }),
  }))
}

fn fixed_timescale(tc: f64, graph: &Graph, node_times: &CoalescentNodeTimes) -> Result<CoalescentTimescale, Report> {
  let lineage_counts =
    compute_lineage_counts(graph, node_times).wrap_err("Failed to compute coalescent lineage counts")?;
  let breakpoints = lineage_counts.breakpoints();
  if breakpoints.is_empty() {
    return make_error!("Cannot report a fixed coalescent Tc: the tree has no node times to span");
  }
  let t_min = breakpoints[0];
  let t_max = breakpoints[breakpoints.len() - 1];
  Ok(CoalescentTimescale {
    report: Some(CoalescentTcReport {
      segment_boundaries: array![t_min, t_max],
      band: None,
      log_likelihood: None,
    }),
    ..CoalescentTimescale::constant(tc)
  })
}

pub fn build_coalescent_output(
  requested: CoalescentMode,
  timescale: &CoalescentTimescale,
  gen_per_year: f64,
  skyline_params: &SkylineParams,
) -> Result<Option<CoalescentOutput>, Report> {
  let Some(mode) = requested.output_mode() else {
    return Ok(None);
  };
  let report = timescale.report.as_ref().ok_or_else(|| {
    make_report!("A coalescent output ({mode:?}) must carry a per-segment report, but none was produced")
  })?;

  let tc_values = timescale.schedule.values().to_vec();
  let boundaries = report.segment_boundaries.to_vec();

  let (n_points, stiffness) = match mode {
    CoalescentOutputMode::Skyline => (Some(skyline_params.n_points), Some(skyline_params.stiffness)),
    CoalescentOutputMode::Fixed | CoalescentOutputMode::Constant => (None, None),
  };
  let confidence_n_std = report.band.as_ref().map(|_| skyline_params.n_std);

  let (lower, upper) = match &report.band {
    Some(band) => (band.lower.to_vec(), band.upper.to_vec()),
    None => (Vec::new(), Vec::new()),
  };
  let band = report.band.as_ref().map(|_| CoalescentBand {
    lower: &lower,
    upper: &upper,
  });

  let output = CoalescentOutput::new(
    CoalescentInputs {
      mode,
      n_points,
      stiffness,
      confidence_n_std,
      gen_per_year,
    },
    &CoalescentSolve {
      segment_boundaries: &boundaries,
      tc_values: &tc_values,
      band,
      log_likelihood: report.log_likelihood,
    },
  )?;
  Ok(Some(output))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoalescentMode {
  Disabled,
  Fixed(f64),
  Constant,
  Skyline,
}

impl CoalescentMode {
  pub fn is_optimized(self) -> bool {
    matches!(self, CoalescentMode::Constant | CoalescentMode::Skyline)
  }

  pub fn output_mode(self) -> Option<CoalescentOutputMode> {
    match self {
      CoalescentMode::Disabled => None,
      CoalescentMode::Fixed(_) => Some(CoalescentOutputMode::Fixed),
      CoalescentMode::Constant => Some(CoalescentOutputMode::Constant),
      CoalescentMode::Skyline => Some(CoalescentOutputMode::Skyline),
    }
  }
}

pub struct CoalescentTimescale {
  pub distribution: Distribution,
  pub schedule: PiecewiseConstantFn,
  pub report: Option<CoalescentTcReport>,
}

impl CoalescentTimescale {
  pub fn constant(tc: f64) -> Self {
    Self {
      distribution: Distribution::constant(tc),
      schedule: PiecewiseConstantFn::new(array![], array![tc]),
      report: None,
    }
  }
}

pub struct CoalescentTcReport {
  pub segment_boundaries: Array1<f64>,
  pub band: Option<CoalescentReportBand>,
  pub log_likelihood: Option<f64>,
}

pub struct CoalescentReportBand {
  pub lower: Array1<f64>,
  pub upper: Array1<f64>,
}
