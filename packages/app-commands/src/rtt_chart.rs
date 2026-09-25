use crate::rtt_chart_render::draw_chart;
use eyre::Report;
use itertools::{Itertools, chain};
use plotters::prelude::*;
use std::path::Path;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::rtt::ClockRegressionResult;
use treetime::progress::ProgressSink;
#[cfg(not(feature = "png"))]
use treetime::progress_warn;

#[cfg(feature = "png")]
use image::{ColorType, DynamicImage, ImageBuffer, ImageEncoder, Rgb, codecs::png::PngEncoder};
#[cfg(feature = "png")]
use treetime_utils::io::file::create_file_or_stdout;

const CHART_SIZE: (u32, u32) = (1200, 800);

pub fn write_clock_regression_chart_svg(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  filepath: impl AsRef<Path>,
) -> Result<(), Report> {
  let svg = SVGBackend::new(filepath.as_ref(), CHART_SIZE).into_drawing_area();
  draw_chart(results, clock_model, &svg)?;
  svg.present()?;
  Ok(())
}

#[cfg(feature = "png")]
pub fn write_clock_regression_chart_png(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  filepath: impl AsRef<Path>,
  _progress: &dyn ProgressSink,
) -> Result<(), Report> {
  let img = write_clock_regression_chart_bitmap(results, clock_model)?;
  let mut f = &mut create_file_or_stdout(filepath)?;
  let encoder = PngEncoder::new(&mut f);
  encoder.write_image(img.as_bytes(), img.width(), img.height(), ColorType::Rgb8.into())?;
  Ok(())
}

#[cfg(not(feature = "png"))]
pub fn write_clock_regression_chart_png(
  _results: &[ClockRegressionResult],
  _clock_model: &ClockModel,
  _filepath: impl AsRef<Path>,
  progress: &dyn ProgressSink,
) -> Result<(), Report> {
  progress_warn!(
    progress,
    "PNG chart output requested but binary was built without the 'png' feature"
  );
  Ok(())
}

#[cfg(feature = "png")]
fn write_clock_regression_chart_bitmap(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
) -> Result<DynamicImage, Report> {
  let (width, height) = CHART_SIZE;
  let mut img: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::new(width, height);
  {
    let bitmap = BitMapBackend::with_buffer(img.as_flat_samples_mut().samples, (width, height)).into_drawing_area();
    draw_chart(results, clock_model, &bitmap)?;
    bitmap.present()?;
  }
  Ok(DynamicImage::ImageRgb8(img))
}

#[allow(
  clippy::as_conversions,
  clippy::unwrap_used,
  reason = "count/index numeric cast is exact for the domain range; unwrap on a value an upstream invariant guarantees is present"
)]
pub fn gather_points(results: &[ClockRegressionResult], clock_model: &ClockModel) -> Result<PointsResult, Report> {
  assert!(!results.is_empty());

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

  let (x_min, x_max) = points.iter().map(|(x, _)| *x).minmax().into_option().unwrap();

  let line_y1 = clock_model.div(x_min as f64) as f32;
  let line_y2 = clock_model.div(x_max as f64) as f32;
  let line = [(x_min, line_y1), (x_max, line_y2)];

  let (y_min, y_max) = points
    .iter()
    .map(|(_, y)| *y)
    .chain([line_y1, line_y2])
    .minmax()
    .into_option()
    .unwrap();

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
