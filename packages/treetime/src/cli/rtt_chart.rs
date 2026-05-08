use crate::cli::rtt_chart_render::draw_chart;
use crate::clock::clock_model::ClockModel;
use crate::clock::rtt::ClockRegressionResult;
use crate::o;
use comfy_table::modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS};
use comfy_table::presets::UTF8_FULL;
use comfy_table::{ContentArrangement, Table};
use crossterm::terminal;
use eyre::Report;
use itertools::{Itertools, chain};
use num_traits::clamp;
use plotters::prelude::*;
use rgb::RGB8;
use std::path::Path;
use textplots::{Chart, ColorPlot, Shape};

#[cfg(feature = "png")]
use image::{ColorType, DynamicImage, ImageBuffer, ImageEncoder, Rgb, codecs::png::PngEncoder};
#[cfg(feature = "png")]
use treetime_utils::io::file::create_file_or_stdout;

pub fn write_clock_regression_chart_svg(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  filepath: impl AsRef<Path>,
) -> Result<(), Report> {
  let svg = SVGBackend::new(filepath.as_ref(), (1200, 800)).into_drawing_area();
  draw_chart(results, clock_model, &svg)?;
  svg.present()?;
  Ok(())
}

#[cfg(feature = "png")]
pub fn write_clock_regression_chart_png(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  filepath: impl AsRef<Path>,
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
) -> Result<(), Report> {
  log::warn!("PNG chart output requested but binary was built without the 'png' feature");
  Ok(())
}

#[cfg(feature = "png")]
pub fn write_clock_regression_chart_bitmap(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
) -> Result<DynamicImage, Report> {
  let (width, height) = (1200, 800);
  let mut img: ImageBuffer<Rgb<u8>, Vec<u8>> = ImageBuffer::new(width, height);
  {
    let bitmap = BitMapBackend::with_buffer(img.as_flat_samples_mut().samples, (width, height)).into_drawing_area();
    draw_chart(results, clock_model, &bitmap)?;
    bitmap.present()?;
  }
  Ok(DynamicImage::ImageRgb8(img))
}

pub fn print_clock_regression_chart(results: &[ClockRegressionResult], clock_model: &ClockModel) -> Result<(), Report> {
  let mut table = Table::new();
  table
    .load_preset(UTF8_FULL)
    .apply_modifier(UTF8_ROUND_CORNERS)
    .apply_modifier(UTF8_SOLID_INNER_BORDERS)
    .set_content_arrangement(ContentArrangement::Dynamic);

  table.add_row([o!("Clock regression"), clock_model.equation_str()]);
  table.add_row([o!("tMRCA"), format!("{:.1}", clock_model.t_mrca())]);
  table.add_row([o!("Rate"), format!("{:.6}", clock_model.clock_rate())]);
  table.add_row([o!("Intercept"), format!("{:.4}", clock_model.intercept())]);
  if let Some(r_val) = clock_model.r_val() {
    table.add_row([o!("R"), format!("{r_val:.4}")]);
    table.add_row([o!("R²"), format!("{:.4}", r_val.powf(2.0))]);
  }
  if let Some(chisq) = clock_model.chisq() {
    table.add_row([o!("χ²"), format!("{chisq:.3e}")]);
  }
  println!("{table}");

  let (width, height) = terminal::size().unwrap_or((120, 40));
  let width = clamp(width, 0, 1024) as u32;
  let height = clamp(height, 0, 1024) as u32;

  let PointsResult {
    norm_points,
    outlier_points,
    x_min,
    x_max,
    y_min,
    y_max,
    ..
  } = gather_points(results, clock_model)?;

  let mut chart = Chart::new_with_y_range(width, height, x_min, x_max, y_min, y_max);

  let norm_points = Shape::Points(&norm_points);
  let chart = chart.linecolorplot(&norm_points, RGB8 { r: 8, g: 232, b: 140 });

  let outlier_points = Shape::Points(&outlier_points);
  let chart = chart.linecolorplot(&outlier_points, RGB8 { r: 255, g: 105, b: 97 });

  let line = Box::new(|date: f32| clock_model.div(date as f64) as f32);
  let line = Shape::Continuous(line);
  let chart = chart.linecolorplot(&line, RGB8 { r: 8, g: 140, b: 232 });

  chart.nice();

  Ok(())
}

#[allow(clippy::field_scoped_visibility_modifiers)]
pub(crate) struct PointsResult {
  pub(crate) norm_points: Vec<(f32, f32)>,
  pub(crate) outlier_points: Vec<(f32, f32)>,
  pub(crate) line: [(f32, f32); 2],
  pub(crate) x_min: f32,
  pub(crate) x_max: f32,
  pub(crate) y_min: f32,
  pub(crate) y_max: f32,
}

pub(crate) fn gather_points(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
) -> Result<PointsResult, Report> {
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
