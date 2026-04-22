use crate::testing::framework::results::TestResult;
use crate::testing::framework::test_case::TestCase;
use crate::testing::plots::utils::expand_range;
use eyre::Report;
use plotters::prelude::*;
use std::iter::once;

pub fn plot_pointwise_error_profiles<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let output_path = format!("{output_dir}/pointwise_error_profiles.svg");
  let root = SVGBackend::new(&output_path, (1200, 900)).into_drawing_area();
  root.fill(&WHITE)?;

  let areas = root.split_evenly((2, 2));
  let pw = &result.metrics.pointwise;
  let x = &result.evaluation_grid;

  {
    let (x_min, x_max) = expand_range(
      x.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
      x.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),
    );
    let (y_min, y_max) = expand_range(0.0, pw.errors.absolute.iter().fold(0.0_f64, |a, &b| a.max(b)));
    let mut chart = ChartBuilder::on(&areas[0])
      .caption("Absolute Error", ("Arial", 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
      x.iter().zip(pw.errors.absolute.iter()).map(|(&x, &y)| (x, y)),
      RED.stroke_width(2),
    ))?;
  }

  {
    let (x_min, x_max) = expand_range(
      x.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
      x.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),
    );
    let max_rel = pw
      .errors
      .relative
      .iter()
      .copied()
      .filter(|v| v.is_finite())
      .fold(0.0_f64, f64::max);
    let (y_min, y_max) = expand_range(0.0, max_rel);
    let mut chart = ChartBuilder::on(&areas[1])
      .caption("Relative Error", ("Arial", 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
      x.iter()
        .zip(pw.errors.relative.iter())
        .filter_map(|(&x, &y)| y.is_finite().then_some((x, y))),
      BLUE.stroke_width(2),
    ))?;
  }

  {
    let (x_min, x_max) = expand_range(
      x.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
      x.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),
    );
    let min_signed = pw.errors.signed.iter().copied().fold(f64::INFINITY, f64::min);
    let max_signed = pw.errors.signed.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let (y_min, y_max) = expand_range(min_signed, max_signed);
    let mut chart = ChartBuilder::on(&areas[2])
      .caption("Signed Error (bias)", ("Arial", 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
      x.iter().zip(pw.errors.signed.iter()).map(|(&x, &y)| (x, y)),
      GREEN.stroke_width(2),
    ))?;
    chart.draw_series(once(PathElement::new(
      vec![(x_min, 0.0), (x_max, 0.0)],
      BLACK.stroke_width(1),
    )))?;
  }

  {
    let (x_min, x_max) = expand_range(
      x.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
      x.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),
    );
    let max_log = pw
      .errors
      .logarithmic
      .iter()
      .copied()
      .filter(|v| v.is_finite())
      .fold(0.0_f64, f64::max);
    let (y_min, y_max) = expand_range(0.0, max_log);
    let mut chart = ChartBuilder::on(&areas[3])
      .caption("Logarithmic Error (tails)", ("Arial", 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
      x.iter()
        .zip(pw.errors.logarithmic.iter())
        .filter_map(|(&x, &y)| y.is_finite().then_some((x, y))),
      MAGENTA.stroke_width(2),
    ))?;
  }

  root.present()?;
  Ok(())
}

pub fn plot_derivative_errors<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let output_path = format!("{output_dir}/derivative_errors.svg");
  let root = SVGBackend::new(&output_path, (1200, 450)).into_drawing_area();
  root.fill(&WHITE)?;

  let areas = root.split_evenly((1, 2));
  let pw = &result.metrics.pointwise;
  let x = &result.evaluation_grid;
  let (x_min, x_max) = expand_range(
    x.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
    x.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),
  );

  {
    let max_d1 = pw.structural.first_derivative.iter().copied().fold(0.0_f64, f64::max);
    let (y_min, y_max) = expand_range(0.0, max_d1);
    let mut chart = ChartBuilder::on(&areas[0])
      .caption("First Derivative Error", ("Arial", 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
      x.iter()
        .zip(pw.structural.first_derivative.iter())
        .map(|(&x, &y)| (x, y)),
      RED.stroke_width(2),
    ))?;
  }

  {
    let max_d2 = pw.structural.second_derivative.iter().copied().fold(0.0_f64, f64::max);
    let (y_min, y_max) = expand_range(0.0, max_d2);
    let mut chart = ChartBuilder::on(&areas[1])
      .caption("Second Derivative Error", ("Arial", 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().draw()?;
    chart.draw_series(LineSeries::new(
      x.iter()
        .zip(pw.structural.second_derivative.iter())
        .map(|(&x, &y)| (x, y)),
      BLUE.stroke_width(2),
    ))?;
  }

  root.present()?;
  Ok(())
}
