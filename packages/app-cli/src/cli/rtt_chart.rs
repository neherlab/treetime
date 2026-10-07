use app_commands::rtt_chart::{PointsResult, gather_points};
use comfy_table::modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS};
use comfy_table::presets::UTF8_FULL;
use comfy_table::{ContentArrangement, Table};
use crossterm::terminal;
use eyre::{Report, WrapErr};
use num_traits::clamp;
use rgb::RGB8;
use std::io::{self, Write};
use textplots::{Chart, ColorPlot, Shape};
use treetime::clock::clock_model::ClockModel;
use treetime::clock::rtt::ClockRegressionResult;
use treetime::o;

const FALLBACK_TERMINAL_SIZE: (u16, u16) = (120, 40);
const TEXT_CHART_MIN_SIZE: (u16, u16) = (32, 3);
const TEXT_CHART_MAX_SIZE: (u16, u16) = (1024, 1024);

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    result_defaulted,
    reason = "output that is not a terminal has no size; the chart uses a fixed default"
  )
)]
pub(crate) fn print_clock_regression_chart(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
) -> Result<(), Report> {
  let terminal_size = terminal::size()
    .ok()
    .filter(|&(width, height)| width > 0 && height > 0)
    .unwrap_or(FALLBACK_TERMINAL_SIZE);
  write_clock_regression_chart_text(&mut io::stderr().lock(), results, clock_model, terminal_size)
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(crate) fn write_clock_regression_chart_text(
  writer: &mut impl Write,
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  (width, height): (u16, u16),
) -> Result<(), Report> {
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
  if let Some((r_val, r_squared)) = clock_model.r_val().zip(clock_model.r_squared()) {
    table.add_row([o!("R"), format!("{r_val:.4}")]);
    table.add_row([o!("R²"), format!("{r_squared:.4}")]);
  }
  if let Some(chisq) = clock_model.chisq() {
    table.add_row([o!("χ²"), format!("{chisq:.3e}")]);
  }
  writeln!(writer, "{table}").wrap_err("When writing the clock model table")?;

  let width = u32::from(clamp(width, TEXT_CHART_MIN_SIZE.0, TEXT_CHART_MAX_SIZE.0));
  let height = u32::from(clamp(height, TEXT_CHART_MIN_SIZE.1, TEXT_CHART_MAX_SIZE.1));

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

  chart.borders();
  chart.axis();
  chart.figures();
  writeln!(writer, "{chart}").wrap_err("When writing the clock regression chart")?;

  Ok(())
}
