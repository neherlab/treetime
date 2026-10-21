use eyre::Report;
use itertools::{Itertools, izip};
use std::iter::once;
use treetime_utils::fmt::string::counted;
use treetime_utils::make_error;

#[derive(Clone, Debug, PartialEq)]
pub enum TcSchedule {
  Piecewise(PiecewiseTc),
  Exponential(ExponentialTc),
}

impl TcSchedule {
  pub fn constant(tc: f64) -> Result<Self, Report> {
    Self::piecewise(vec![], vec![tc])
  }

  pub fn piecewise(breaks: Vec<f64>, values: Vec<f64>) -> Result<Self, Report> {
    if values.len() != breaks.len() + 1 {
      return make_error!(
        "A piecewise Tc schedule with {} needs {}, got {}",
        counted(breaks.len(), "break"),
        counted(breaks.len() + 1, "value"),
        values.len()
      );
    }
    if let Some(value) = breaks.iter().find(|value| !value.is_finite()) {
      return make_error!("Breaks of a piecewise Tc schedule must be finite, got {value}");
    }
    if !breaks.iter().tuple_windows().all(|(left, right)| left < right) {
      return make_error!("Breaks of a piecewise Tc schedule must be strictly increasing, got {breaks:?}");
    }
    for &tc in &values {
      validate_tc(tc)?;
    }
    Ok(Self::Piecewise(PiecewiseTc { breaks, values }))
  }

  pub fn exponential(tc0: f64, scale: f64) -> Result<Self, Report> {
    validate_tc(tc0)?;
    if !(scale.is_finite() && scale > 0.0) {
      return make_error!("Time scale of an exponential Tc schedule must be finite and positive, got {scale}");
    }
    Ok(Self::Exponential(ExponentialTc { tc0, scale }))
  }

  pub fn tc(&self, tau: f64) -> f64 {
    match self {
      Self::Piecewise(piecewise) => piecewise.values[piecewise.breaks.partition_point(|&value| value <= tau)],
      Self::Exponential(exponential) => exponential.tc0 * (-tau / exponential.scale).exp(),
    }
  }

  pub fn cumulative_hazard(&self, tau: f64) -> f64 {
    match self {
      Self::Piecewise(piecewise) => piecewise
        .segments()
        .take_while(|segment| segment.start < tau)
        .map(|segment| (tau.min(segment.end) - segment.start) / segment.tc)
        .sum(),
      Self::Exponential(exponential) => exponential.scale / exponential.tc0 * (tau / exponential.scale).exp_m1(),
    }
  }

  pub fn inverse_cumulative_hazard(&self, hazard: f64) -> f64 {
    match self {
      Self::Piecewise(piecewise) => {
        let mut hazard_before = 0.0;
        for segment in piecewise.segments() {
          let segment_hazard = (segment.end - segment.start) / segment.tc;
          if hazard <= hazard_before + segment_hazard {
            return segment.start + (hazard - hazard_before) * segment.tc;
          }
          hazard_before += segment_hazard;
        }
        f64::INFINITY
      },
      Self::Exponential(exponential) => exponential.scale * (hazard * exponential.tc0 / exponential.scale).ln_1p(),
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct PiecewiseTc {
  breaks: Vec<f64>,
  values: Vec<f64>,
}

impl PiecewiseTc {
  fn segments(&self) -> impl Iterator<Item = Segment> + '_ {
    let starts = once(f64::NEG_INFINITY).chain(self.breaks.iter().copied());
    let ends = self.breaks.iter().copied().chain([f64::INFINITY]);
    izip!(starts, ends, self.values.iter().copied())
      .filter(|&(_, end, _)| end > 0.0)
      .map(|(start, end, tc)| Segment {
        start: start.max(0.0),
        end,
        tc,
      })
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ExponentialTc {
  tc0: f64,
  scale: f64,
}

fn validate_tc(tc: f64) -> Result<(), Report> {
  if !(tc.is_finite() && tc > 0.0) {
    return make_error!("Coalescent time scale Tc must be finite and positive, got {tc}");
  }
  Ok(())
}

struct Segment {
  start: f64,
  end: f64,
  tc: f64,
}
