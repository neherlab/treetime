use crate::distribution::reference::convolution_test::framework::results::TestResult;
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::plots::utils::{combined_range, expand_range, tolerance_label};
use eyre::Report;
use itertools::izip;
use plotters::prelude::*;

pub fn plot_absolute_error<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let mut error_points = Vec::with_capacity(result.evaluation_grid.len());
  let mut max_error = 0.0_f64;
  for (&x, &actual, &expected) in izip!(
    result.evaluation_grid.iter(),
    result.actual_values.iter(),
    result.expected_values.iter()
  ) {
    let value = (actual - expected).abs();
    if value > max_error {
      max_error = value;
    }
    error_points.push((x, value));
  }
  let (x_raw_min, x_raw_max) = combined_range(&[&result.evaluation_grid]);
  let (x_min, x_max) = expand_range(x_raw_min, x_raw_max);
  let (mut y_min, mut y_max) = expand_range(0.0, max_error);
  if y_min < 0.0 {
    y_min = 0.0;
  }
  if y_max <= y_min {
    y_max = y_min + 1.0;
  }
  let output_path = format!("{output_dir}/absolute_error.svg");
  let root = SVGBackend::new(&output_path, (960, 480)).into_drawing_area();
  root.fill(&WHITE)?;
  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!("{} | {} | Absolute Error", result.test_case.name(), result.algorithm),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(40)
    .y_label_area_size(60)
    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
  chart.configure_mesh().y_desc("|actual - expected|").draw()?;
  chart
    .draw_series(LineSeries::new(error_points.iter().copied(), RED.stroke_width(2)))?
    .label("Absolute error")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], RED.stroke_width(2)));
  chart
    .configure_series_labels()
    .border_style(BLACK)
    .background_style(WHITE.mix(0.8))
    .draw()?;
  root.present()?;
  Ok(())
}

pub fn plot_tolerance_metrics<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let abs = result.metrics.aggregate.domain_agreement.abs_tolerance_fractions;
  let rel = result.metrics.aggregate.domain_agreement.rel_tolerance_fractions;
  let output_path = format!("{output_dir}/tolerance_metrics.svg");
  let root = SVGBackend::new(&output_path, (800, 600)).into_drawing_area();
  root.fill(&WHITE)?;
  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!(
        "{} | {} | Tolerance Fractions",
        result.test_case.name(),
        result.algorithm
      ),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(50)
    .y_label_area_size(60)
    .build_cartesian_2d(-0.2..3.0, 0.0..1.05)?;
  chart
    .configure_mesh()
    .x_labels(3)
    .y_labels(6)
    .x_desc("Tolerance level")
    .y_desc("Fraction within threshold")
    .x_label_formatter(&|value| tolerance_label(*value))
    .y_label_formatter(&|value| format!("{value:.2}"))
    .draw()?;
  chart
    .draw_series((0..3).map(|i| {
      let x0 = i as f64 + 0.05;
      let x1 = x0 + 0.35;
      Rectangle::new([(x0, 0.0), (x1, abs[i])], BLUE.mix(0.6).filled())
    }))?
    .label("Absolute")
    .legend(|(x, y)| Rectangle::new([(x - 6, y - 6), (x + 6, y + 6)], BLUE.mix(0.6).filled()));
  chart
    .draw_series((0..3).map(|i| {
      let x0 = i as f64 + 0.5;
      let x1 = x0 + 0.35;
      Rectangle::new([(x0, 0.0), (x1, rel[i])], RED.mix(0.6).filled())
    }))?
    .label("Relative")
    .legend(|(x, y)| Rectangle::new([(x - 6, y - 6), (x + 6, y + 6)], RED.mix(0.6).filled()));
  chart
    .configure_series_labels()
    .border_style(BLACK)
    .background_style(WHITE.mix(0.8))
    .draw()?;
  root.present()?;
  Ok(())
}

pub fn plot_error_histogram<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let output_path = format!("{output_dir}/error_histogram.svg");
  let root = SVGBackend::new(&output_path, (960, 640)).into_drawing_area();
  root.fill(&WHITE)?;

  let dist = &result.metrics.distribution;
  let hist = &dist.histograms.abs_error_histogram;

  if hist.bin_centers.is_empty() {
    return Ok(());
  }

  let x_min = hist.bin_centers.iter().fold(f64::INFINITY, |a, &b| a.min(b));
  let x_max = hist.bin_centers.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
  let y_max = *hist.bin_counts.iter().max().unwrap_or(&0) as f64;

  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!(
        "{} | {} | Error Distribution",
        result.test_case.name(),
        result.algorithm
      ),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(40)
    .y_label_area_size(60)
    .build_cartesian_2d(x_min..x_max, 0.0..y_max * 1.1)?;

  chart
    .configure_mesh()
    .x_desc("log10(Absolute Error)")
    .y_desc("Count")
    .draw()?;

  let bar_width = if hist.bin_centers.len() > 1 {
    (hist.bin_centers[1] - hist.bin_centers[0]) * 0.8
  } else {
    0.5
  };

  chart.draw_series(
    hist
      .bin_centers
      .iter()
      .zip(hist.bin_counts.iter())
      .map(|(&center, &count)| {
        Rectangle::new(
          [
            (center - bar_width / 2.0, 0.0),
            (center + bar_width / 2.0, count as f64),
          ],
          BLUE.mix(0.6).filled(),
        )
      }),
  )?;

  root.present()?;
  Ok(())
}
