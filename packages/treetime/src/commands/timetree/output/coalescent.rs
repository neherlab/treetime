use crate::coalescent::population_size::effective_population_size;
use eyre::{Report, WrapErr};
use serde::{Deserialize, Serialize};
use std::path::Path;
use treetime_io::csv::{CsvStructFileWriter, CsvStructWriter};
use treetime_utils::io::json::{JsonPretty, json_write_file, json_write_str};
use treetime_utils::{make_error, make_report};

/// Serialize the coalescent output as one rich JSON document (`inputs` + `outputs.segments`).
pub fn coalescent_json_str(output: &CoalescentOutput) -> Result<String, Report> {
  json_write_str(output, JsonPretty(true))
}

/// Serialize the coalescent segments as flat delimited rows (header + one row per segment).
///
/// The header is derived from the serde field names of [`CoalescentSegmentRow`], so the dotted
/// names (`segment.start`, `T_c.value`, ...) become the column titles. `delimiter` selects CSV
/// (`b','`) or TSV (`b'\t'`).
pub fn coalescent_delimited_str(output: &CoalescentOutput, delimiter: u8) -> Result<String, Report> {
  let mut writer = CsvStructWriter::new(Vec::<u8>::new(), delimiter)?;
  for row in output.rows() {
    writer.write(&row)?;
  }
  let bytes = writer
    .writer
    .into_inner()
    .map_err(|err| make_report!("Coalescent delimited serialization failed to flush: {err}"))?;
  String::from_utf8(bytes).wrap_err("Coalescent delimited output is not valid UTF-8")
}

/// Write the rich coalescent document to a JSON file.
pub fn write_coalescent_json(output: &CoalescentOutput, path: impl AsRef<Path>) -> Result<(), Report> {
  json_write_file(path, output, JsonPretty(true))
}

/// Write the flat coalescent segments to a delimited file (CSV with `b','`, TSV with `b'\t'`).
pub fn write_coalescent_delimited(
  output: &CoalescentOutput,
  path: impl AsRef<Path>,
  delimiter: u8,
) -> Result<(), Report> {
  let mut writer = CsvStructFileWriter::new(path, delimiter)?;
  for row in output.rows() {
    writer.write(&row)?;
  }
  Ok(())
}

/// The rich coalescent output document: the run parameters that produced the result and the
/// per-segment time scales.
///
/// This is the single source of truth for both serializations. JSON emits the whole document with
/// nested `segment`/`T_c`/`N_e` objects; the delimited formats emit only [`Self::rows`], the flat
/// per-segment projection.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentOutput {
  pub inputs: CoalescentInputs,
  pub outputs: CoalescentOutputs,
}

impl CoalescentOutput {
  /// Builds the document from a coalescent solve.
  ///
  /// Consumes the optimized per-segment `T_c` values and the P1 confidence band, and derives the
  /// effective population size `N_e = T_c * gen_per_year` (P2) for the point and, when present, for
  /// each band bound. A fixed (user-supplied) `T_c` carries no band, so both `T_c` and `N_e` report
  /// a value only.
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
      outputs: CoalescentOutputs { segments },
    })
  }

  /// Projects the rich segments to flat rows: the documented rich -> flat mapping.
  ///
  /// Each nested [`CoalescentSegment`] becomes one [`CoalescentSegmentRow`]: the 0-based JSON
  /// `index` becomes 1-based, and the nested `segment`/`T_c`/`N_e` objects flatten to dotted
  /// scalar columns. The projection is total, so `rows().len() == outputs.segments.len()`.
  pub fn rows(&self) -> Vec<CoalescentSegmentRow> {
    self
      .outputs
      .segments
      .iter()
      .map(CoalescentSegmentRow::from_segment)
      .collect()
  }
}

/// The run parameters that produced a coalescent result.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentInputs {
  pub mode: CoalescentMode,
  /// Generations per year used to map `T_c` to `N_e`.
  pub gen_per_year: f64,
  /// Confidence level, in standard deviations, of the `T_c` band. `None` for a fixed `T_c`, which
  /// carries no band.
  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence_n_std: Option<f64>,
}

/// How the reported `T_c` was produced. The disabled mode reports no coalescent, so it has no
/// variant here.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CoalescentMode {
  /// Fixed, user-supplied `T_c` (one segment, no band).
  Fixed,
  /// Optimized constant `T_c` (one segment).
  Constant,
  /// Optimized piecewise-constant skyline `T_c(t)`.
  Skyline,
}

/// The per-segment coalescent time scales.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentOutputs {
  pub segments: Vec<CoalescentSegment>,
}

/// One coalescent segment in rich (nested) form, as emitted to JSON.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentSegment {
  /// 0-based segment index.
  pub index: usize,
  pub segment: SegmentInterval,
  #[serde(rename = "T_c")]
  pub tc: Estimate,
  #[serde(rename = "N_e")]
  pub ne: Estimate,
}

/// One coalescent segment in flat form, as emitted to CSV/TSV. Serde field renames produce the
/// dotted column titles that mirror the nested JSON keys.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct CoalescentSegmentRow {
  /// 1-based segment index.
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

impl CoalescentSegmentRow {
  fn from_segment(segment: &CoalescentSegment) -> Self {
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

/// A half-open numeric-date span `[start, end)` for one coalescent segment.
///
/// Local stand-in for the shared numeric-date interval type. The wider unification onto
/// `DateRangeNumeric` is deferred (see `kb/issues/N-datetime-date-and-range-representation-inconsistent.md`).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct SegmentInterval {
  pub start: f64,
  pub end: f64,
}

/// A point estimate with an optional confidence band. `lower`/`upper` are absent only for a fixed
/// `T_c` (and the `N_e` derived from it), which carries no band.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Estimate {
  pub value: f64,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub lower: Option<f64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub upper: Option<f64>,
}

impl Estimate {
  /// A bare point estimate with no band.
  pub fn point(value: f64) -> Self {
    Self {
      value,
      lower: None,
      upper: None,
    }
  }

  /// A point estimate with a confidence band.
  pub fn with_band(value: f64, lower: f64, upper: f64) -> Self {
    Self {
      value,
      lower: Some(lower),
      upper: Some(upper),
    }
  }
}

/// A coalescent solve to serialize: the optimized per-segment `T_c` and its segment boundaries,
/// plus the optional P1 confidence band.
pub struct CoalescentSolve<'a> {
  /// Segment boundaries in numeric date (length `tc_values.len() + 1`, ascending).
  pub segment_boundaries: &'a [f64],
  /// Optimized `T_c` per segment.
  pub tc_values: &'a [f64],
  /// Per-segment `T_c` confidence band, or `None` for a fixed `T_c`.
  pub band: Option<CoalescentBand<'a>>,
}

/// Per-segment `T_c` confidence band bounds. Both bounds are present together, so a partial band is
/// unrepresentable.
pub struct CoalescentBand<'a> {
  pub lower: &'a [f64],
  pub upper: &'a [f64],
}
