use crate::distribution::reference::convolution_test::framework::results::TestResult;
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::plots::utils::{combined_range, expand_range};
use eyre::Report;
use plotters::prelude::*;

pub fn plot_functions_and_convolution<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let output_path = format!("{output_dir}/functions_convolution.svg");
  let root = SVGBackend::new(&output_path, (960, 640)).into_drawing_area();
  root.fill(&WHITE)?;
  let (x_raw_min, x_raw_max) = combined_range(&[&result.f_x_values, &result.g_x_values, &result.evaluation_grid]);
  let (x_min, x_max) = expand_range(x_raw_min, x_raw_max);
  let (y_raw_min, y_raw_max) = combined_range(&[
    &result.f_y_values,
    &result.g_y_values,
    &result.actual_values,
    &result.expected_values,
  ]);
  let (y_min, y_max) = expand_range(y_raw_min, y_raw_max);
  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!(
        "{} | {} | Input and Convolution",
        result.test_case.name(),
        result.algorithm
      ),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(40)
    .y_label_area_size(60)
    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
  chart.configure_mesh().draw()?;
  chart
    .draw_series(LineSeries::new(
      result
        .f_x_values
        .iter()
        .zip(result.f_y_values.iter())
        .map(|(&x, &y)| (x, y)),
      BLUE.stroke_width(2),
    ))?
    .label("f(x)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], BLUE.stroke_width(2)));
  chart
    .draw_series(LineSeries::new(
      result
        .g_x_values
        .iter()
        .zip(result.g_y_values.iter())
        .map(|(&x, &y)| (x, y)),
      GREEN.stroke_width(2),
    ))?
    .label("g(x)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], GREEN.stroke_width(2)));
  chart
    .draw_series(LineSeries::new(
      result
        .evaluation_grid
        .iter()
        .zip(result.actual_values.iter())
        .map(|(&x, &y)| (x, y)),
      MAGENTA.stroke_width(3),
    ))?
    .label(format!("Actual ({})", result.algorithm))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], MAGENTA.stroke_width(3)));
  chart
    .draw_series(LineSeries::new(
      result
        .evaluation_grid
        .iter()
        .zip(result.expected_values.iter())
        .map(|(&x, &y)| (x, y)),
      RED.stroke_width(3),
    ))?
    .label("Expected")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], RED.stroke_width(3)));
  chart
    .configure_series_labels()
    .border_style(BLACK)
    .background_style(WHITE.mix(0.8))
    .draw()?;
  root.present()?;
  Ok(())
}
