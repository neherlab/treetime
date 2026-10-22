use crate::rtt_chart_render::draw_chart;
use eyre::{Report, WrapErr};
use itertools::{Itertools, chain};
use plotters::prelude::*;
use std::io::Write;
use std::path::Path;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::rtt::ClockRegressionResult;
use treetime_charts::convert::ChartFormat;
use treetime_charts::font::register_chart_fonts;
use treetime_utils::io::file::write_file_with;
use treetime_utils::make_error;

const CHART_SIZE: (u32, u32) = (1200, 800);

pub fn write_clock_regression_charts(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  outputs: &[(ChartFormat, &Path)],
) -> Result<(), Report> {
  if outputs.is_empty() {
    return Ok(());
  }
  let svg = clock_regression_chart_svg(results, clock_model)?;
  outputs.iter().try_for_each(|&(format, filepath)| {
    let bytes = format.encode(&svg)?;
    write_file_with(filepath, |writer| {
      writer
        .write_all(&bytes)
        .wrap_err_with(|| format!("When writing the {format:?} chart"))
    })
  })
}

pub(crate) fn clock_regression_chart_svg(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
) -> Result<String, Report> {
  register_chart_fonts()?;
  let mut text = String::new();
  {
    let svg = SVGBackend::with_string(&mut text, CHART_SIZE).into_drawing_area();
    draw_chart(results, clock_model, &svg)?;
    svg.present()?;
  }
  Ok(text)
}

#[expect(
  clippy::as_conversions,
  reason = "the chart library draws f32 coordinates, so dates and divergences narrow from f64 to f32"
)]
pub fn gather_points(results: &[ClockRegressionResult], clock_model: &ClockModel) -> Result<PointsResult, Report> {
  let (outliers, norms): (Vec<_>, Vec<_>) = results.iter().partition(|result| result.is_outlier);

  let norm_points = norms
    .into_iter()
    .filter_map(|result| result.date.map(|date| (date as f32, result.div as f32)))
    .collect_vec();

  let outlier_points = outliers
    .into_iter()
    .filter_map(|result| result.date.map(|date| (date as f32, result.div as f32)))
    .collect_vec();

  let points = chain!(&norm_points, &outlier_points).copied().collect_vec();

  let Some((x_min, x_max)) = points.iter().map(|(x, _)| *x).minmax().into_option() else {
    return make_error!("the root-to-tip chart needs at least one sample with a date");
  };

  let line_y1 = clock_model.div(f64::from(x_min)) as f32;
  let line_y2 = clock_model.div(f64::from(x_max)) as f32;
  let line = [(x_min, line_y1), (x_max, line_y2)];

  let (y_min, y_max) = points
    .iter()
    .map(|(_, y)| *y)
    .chain([line_y1, line_y2])
    .minmax()
    .into_option()
    .unwrap_or((line_y1, line_y2));

  Ok(PointsResult {
    norm_points,
    outlier_points,
    line,
    x_min,
    x_max,
    y_min,
    y_max,
  })
}

pub struct PointsResult {
  pub norm_points: Vec<(f32, f32)>,
  pub outlier_points: Vec<(f32, f32)>,
  pub line: [(f32, f32); 2],
  pub x_min: f32,
  pub x_max: f32,
  pub y_min: f32,
  pub y_max: f32,
}
