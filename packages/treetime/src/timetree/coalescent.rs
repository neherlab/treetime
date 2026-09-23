use crate::coalescent::population_size::effective_population_size;
use eyre::Report;
use serde::{Deserialize, Serialize};
use treetime_utils::make_error;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentOutput {
  pub inputs: CoalescentInputs,
  pub outputs: CoalescentOutputs,
}

impl CoalescentOutput {
  pub fn new(inputs: CoalescentInputs, solve: &CoalescentSolve) -> Result<Self, Report> {
    let n = solve.tc_values.len();
    if solve.segment_boundaries.len() != n + 1 {
      return make_error!(
        "Coalescent output expects {} segment boundaries for {n} segment(s), got {}",
        n + 1,
        solve.segment_boundaries.len()
      );
    }
    if let Some(band) = &solve.band {
      if band.lower.len() != n || band.upper.len() != n {
        return make_error!(
          "Coalescent band bounds must have one entry per segment ({n}), got lower={}, upper={}",
          band.lower.len(),
          band.upper.len()
        );
      }
    }

    let gen_per_year = inputs.gen_per_year;
    let segments = (0..n)
      .map(|i| {
        let tc = match &solve.band {
          Some(band) => Estimate::with_band(solve.tc_values[i], band.lower[i], band.upper[i]),
          None => Estimate::point(solve.tc_values[i]),
        };
        let ne = Estimate {
          value: effective_population_size(tc.value, gen_per_year),
          lower: tc.lower.map(|lower| effective_population_size(lower, gen_per_year)),
          upper: tc.upper.map(|upper| effective_population_size(upper, gen_per_year)),
        };
        CoalescentSegment {
          index: i,
          segment: SegmentInterval {
            start: solve.segment_boundaries[i],
            end: solve.segment_boundaries[i + 1],
          },
          tc,
          ne,
        }
      })
      .collect();

    Ok(Self {
      inputs,
      outputs: CoalescentOutputs {
        log_likelihood: solve.log_likelihood,
        segments,
      },
    })
  }

  pub fn rows(&self) -> Vec<CoalescentSegmentRow> {
    self.outputs.segments.iter().map(CoalescentSegmentRow::from).collect()
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentInputs {
  pub mode: CoalescentOutputMode,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub n_points: Option<usize>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub stiffness: Option<f64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence_n_std: Option<f64>,
  pub gen_per_year: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum CoalescentOutputMode {
  Fixed,
  Constant,
  Skyline,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentOutputs {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub log_likelihood: Option<f64>,
  pub segments: Vec<CoalescentSegment>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentSegmentRow {
  pub index: usize,
  #[serde(rename = "segment.start")]
  pub segment_start: f64,
  #[serde(rename = "segment.end")]
  pub segment_end: f64,
  #[serde(rename = "T_c.value")]
  pub tc_value: f64,
  #[serde(rename = "T_c.lower")]
  pub tc_lower: Option<f64>,
  #[serde(rename = "T_c.upper")]
  pub tc_upper: Option<f64>,
  #[serde(rename = "N_e.value")]
  pub ne_value: f64,
  #[serde(rename = "N_e.lower")]
  pub ne_lower: Option<f64>,
  #[serde(rename = "N_e.upper")]
  pub ne_upper: Option<f64>,
}

impl From<&CoalescentSegment> for CoalescentSegmentRow {
  fn from(segment: &CoalescentSegment) -> Self {
    Self {
      index: segment.index + 1,
      segment_start: segment.segment.start,
      segment_end: segment.segment.end,
      tc_value: segment.tc.value,
      tc_lower: segment.tc.lower,
      tc_upper: segment.tc.upper,
      ne_value: segment.ne.value,
      ne_lower: segment.ne.lower,
      ne_upper: segment.ne.upper,
    }
  }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentSegment {
  pub index: usize,
  pub segment: SegmentInterval,
  #[serde(rename = "T_c")]
  pub tc: Estimate,
  #[serde(rename = "N_e")]
  pub ne: Estimate,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SegmentInterval {
  pub start: f64,
  pub end: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Estimate {
  value: f64,
  #[serde(skip_serializing_if = "Option::is_none")]
  lower: Option<f64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  upper: Option<f64>,
}

impl Estimate {
  pub fn point(value: f64) -> Self {
    Self {
      value,
      lower: None,
      upper: None,
    }
  }

  pub fn with_band(value: f64, lower: f64, upper: f64) -> Self {
    Self {
      value,
      lower: Some(lower),
      upper: Some(upper),
    }
  }
}

pub struct CoalescentSolve<'a> {
  pub segment_boundaries: &'a [f64],
  pub tc_values: &'a [f64],
  pub band: Option<CoalescentBand<'a>>,
  pub log_likelihood: Option<f64>,
}

pub struct CoalescentBand<'a> {
  pub lower: &'a [f64],
  pub upper: &'a [f64],
}
