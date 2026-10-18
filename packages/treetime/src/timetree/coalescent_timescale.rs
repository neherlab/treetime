use crate::coalescent::mean_field::LineageExpectation;
use crate::coalescent::skyline::{SkylineParams, expected_log_likelihood, optimize_skyline};
use crate::progress::LogSink;
use crate::timetree::coalescent::{
  CoalescentBand, CoalescentInputs, CoalescentOutput, CoalescentOutputMode, CoalescentSolve,
};
use eyre::{Report, WrapErr};
use ndarray::{Array1, array};
use std::sync::Arc;
use treetime_distribution::Distribution;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_primitives::LogLh;
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

/// The coalescent time scale of `mode` for the lineages of `expectation`: the fixed Tc, or the constant or skyline Tc
/// fitted to them. Without a coalescent model, the constant fit still gives the merger rates of polytomy resolution.
pub(crate) fn coalescent_timescale(
  mode: CoalescentMode,
  skyline_params: &SkylineParams,
  expectation: &LineageExpectation,
  log: &dyn LogSink,
) -> Result<CoalescentTimescale, Report> {
  let n_points = match mode {
    CoalescentMode::Fixed(tc) => return Ok(fixed_timescale(tc, expectation)),
    CoalescentMode::Disabled | CoalescentMode::Constant => 1,
    CoalescentMode::Skyline => skyline_params.n_points,
  };
  let result = optimize_skyline(
    &SkylineParams {
      n_points,
      ..skyline_params.clone()
    },
    expectation,
    log,
  )
  .wrap_err("Failed to estimate the coalescent timescale")?;
  Ok(CoalescentTimescale {
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
  })
}

fn fixed_timescale(tc: f64, expectation: &LineageExpectation) -> CoalescentTimescale {
  let (t_min, t_max) = expectation.span();
  CoalescentTimescale {
    report: Some(CoalescentTcReport {
      segment_boundaries: array![t_min, t_max],
      band: None,
      log_likelihood: None,
    }),
    ..CoalescentTimescale::constant(tc)
  }
}

pub(crate) fn build_coalescent_output(
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

pub(crate) struct CoalescentTimescale {
  pub distribution: Distribution,
  pub schedule: PiecewiseConstantFn,
  pub report: Option<CoalescentTcReport>,
}

impl CoalescentTimescale {
  /// Expected complete-data log-likelihood of the node times behind `expectation` under this time scale (see
  /// `fn expected_log_likelihood`), on the segments of the report or, without one, the span of the expectation
  pub(crate) fn expected_log_likelihood(&self, expectation: &LineageExpectation) -> LogLh {
    let boundaries = self.report.as_ref().map_or_else(
      || {
        let (t_min, t_max) = expectation.span();
        vec![t_min, t_max]
      },
      |report| report.segment_boundaries.to_vec(),
    );
    let (i_seg, m_seg) = expectation.segment_statistics(&boundaries);
    expected_log_likelihood(&i_seg, &m_seg, self.schedule.values().as_slice().unwrap_or_default())
  }

  pub(crate) fn constant(tc: f64) -> Self {
    Self {
      distribution: Distribution::constant(tc),
      schedule: PiecewiseConstantFn::new(array![], array![tc]),
      report: None,
    }
  }
}

pub(crate) struct CoalescentTcReport {
  pub segment_boundaries: Array1<f64>,
  pub band: Option<CoalescentReportBand>,
  pub log_likelihood: Option<f64>,
}

pub(crate) struct CoalescentReportBand {
  pub lower: Array1<f64>,
  pub upper: Array1<f64>,
}

pub(crate) struct CoalescentSetup {
  pub mode: CoalescentMode,
  pub skyline_params: SkylineParams,
  /// Expected lineage counts from node posteriors, for the mean-field prior on node times and the merger rates of
  /// polytomy resolution
  pub lineage_expectation: Arc<LineageExpectation>,
}

impl CoalescentSetup {
  pub(crate) fn prior_wanted(&self) -> bool {
    self.mode.wants_prior()
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoalescentMode {
  Disabled,
  Fixed(f64),
  Constant,
  Skyline,
}

impl CoalescentMode {
  pub(crate) fn wants_prior(self) -> bool {
    self != CoalescentMode::Disabled
  }

  pub(crate) fn is_optimized(self) -> bool {
    matches!(self, CoalescentMode::Constant | CoalescentMode::Skyline)
  }

  pub(crate) fn output_mode(self) -> Option<CoalescentOutputMode> {
    match self {
      CoalescentMode::Disabled => None,
      CoalescentMode::Fixed(_) => Some(CoalescentOutputMode::Fixed),
      CoalescentMode::Constant => Some(CoalescentOutputMode::Constant),
      CoalescentMode::Skyline => Some(CoalescentOutputMode::Skyline),
    }
  }
}
