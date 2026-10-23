use app_commands::rtt_chart::data::{RttChartData, RttPoint};
use app_commands::rtt_chart::frame::RttFrame;
use app_commands::rtt_chart::summary::RttSummary;
use comfy_table::modifiers::{UTF8_ROUND_CORNERS, UTF8_SOLID_INNER_BORDERS};
use comfy_table::presets::UTF8_FULL;
use comfy_table::{ContentArrangement, Table};
use crossterm::terminal;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use num_traits::clamp;
use rgb::RGB8;
use std::io::{self, Write};
use std::iter::once;
use textplots::{Chart, ColorPlot, Shape};
use treetime::clock::clock_model::ClockModel;
use treetime::clock::rtt::ClockRegressionResult;
use treetime::o;

const FALLBACK_TERMINAL_SIZE: (u16, u16) = (120, 40);
const TEXT_CHART_MIN_SIZE: (u16, u16) = (32, 3);
const TEXT_CHART_MAX_SIZE: (u16, u16) = (1024, 1024);
const SAMPLE_COLOR: RGB8 = RGB8 { r: 143, g: 176, b: 236 };
const EXCLUDED_COLOR: RGB8 = RGB8 { r: 240, g: 135, b: 124 };
const FIT_COLOR: RGB8 = RGB8 { r: 92, g: 195, b: 165 };

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

#[expect(
  clippy::as_conversions,
  reason = "the chart library draws f32 coordinates, so dates and divergences narrow from f64 to f32"
)]
pub(crate) fn write_clock_regression_chart_text(
  writer: &mut impl Write,
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  (width, height): (u16, u16),
) -> Result<(), Report> {
  let data = RttChartData::new(results)?;
  let frame = RttFrame::new(&data, clock_model);
  let summary = RttSummary::new(&data, clock_model);

  let mut table = Table::new();
  table
    .load_preset(UTF8_FULL)
    .apply_modifier(UTF8_ROUND_CORNERS)
    .apply_modifier(UTF8_SOLID_INNER_BORDERS)
    .set_content_arrangement(ContentArrangement::Dynamic);
  for statistic in summary.statistics() {
    table.add_row([statistic.label.clone(), statistic.value_text()]);
  }
  let samples = once(&summary.samples).chain(&summary.notes).join(". ");
  table.add_row([o!("Samples"), samples]);
  writeln!(writer, "{table}").wrap_err("When writing the clock model table")?;

  let width = u32::from(clamp(width, TEXT_CHART_MIN_SIZE.0, TEXT_CHART_MAX_SIZE.0));
  let height = u32::from(clamp(height, TEXT_CHART_MIN_SIZE.1, TEXT_CHART_MAX_SIZE.1));
  let narrow = |points: &[RttPoint]| {
    points
      .iter()
      .map(|point| (point.date as f32, point.div as f32))
      .collect_vec()
  };
  let (fitted, excluded) = (narrow(&data.fitted), narrow(&data.excluded));
  let ((x_min, x_max), (y_min, y_max)) = (frame.dates, frame.divergences);

  let mut chart = Chart::new_with_y_range(width, height, x_min as f32, x_max as f32, y_min as f32, y_max as f32);
  let fitted = Shape::Points(&fitted);
  let chart = chart.linecolorplot(&fitted, SAMPLE_COLOR);
  let excluded = Shape::Points(&excluded);
  let chart = chart.linecolorplot(&excluded, EXCLUDED_COLOR);
  let line = Box::new(|date: f32| clock_model.div(f64::from(date)) as f32);
  let line = Shape::Continuous(line);
  let chart = chart.linecolorplot(&line, FIT_COLOR);

  chart.borders();
  chart.axis();
  chart.figures();
  writeln!(writer, "{chart}").wrap_err("When writing the clock regression chart")?;

  Ok(())
}
