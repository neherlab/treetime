use crate::testing::framework::results::TestResult;
use crate::testing::framework::test_case::TestCase;
use crate::testing::plots::utils::{LABEL_FONT, expand_range};
use eyre::Report;
use plotters::coord::Shift;
use plotters::prelude::*;
use std::iter::once;
use treetime_charts::font::CHART_FONT_FAMILY;

pub(crate) fn plot_spatial_profiles<T>(
  result: &TestResult<T>,
  root: &DrawingArea<SVGBackend<'_>, Shift>,
) -> Result<(), Report>
where
  T: TestCase,
{
  root.fill(&WHITE)?;

  let areas = root.split_evenly((2, 2));
  let sp = &result.metrics.spatial;
  let x = &result.evaluation_grid;
  let (x_min, x_max) = expand_range(
    x.iter().fold(f64::INFINITY, |a, &b| a.min(b)),
    x.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b)),
  );

  {
    let max_rms = sp.windowed.sliding_rms.iter().copied().fold(0.0_f64, f64::max);
    let (y_min, y_max) = expand_range(0.0, max_rms);
    let mut chart = ChartBuilder::on(&areas[0])
      .caption("Sliding Window RMS Error", (CHART_FONT_FAMILY, 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().label_style(LABEL_FONT).draw()?;
    chart.draw_series(LineSeries::new(
      x.iter().zip(sp.windowed.sliding_rms.iter()).map(|(&x, &y)| (x, y)),
      RED.stroke_width(2),
    ))?;
  }

  {
    let max_max = sp.windowed.sliding_max.iter().copied().fold(0.0_f64, f64::max);
    let (y_min, y_max) = expand_range(0.0, max_max);
    let mut chart = ChartBuilder::on(&areas[1])
      .caption("Sliding Window Max Error", (CHART_FONT_FAMILY, 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().label_style(LABEL_FONT).draw()?;
    chart.draw_series(LineSeries::new(
      x.iter().zip(sp.windowed.sliding_max.iter()).map(|(&x, &y)| (x, y)),
      BLUE.stroke_width(2),
    ))?;
  }

  {
    let min_cum = sp
      .cumulative
      .cumulative_error
      .iter()
      .copied()
      .fold(f64::INFINITY, f64::min);
    let max_cum = sp
      .cumulative
      .cumulative_error
      .iter()
      .copied()
      .fold(f64::NEG_INFINITY, f64::max);
    let (y_min, y_max) = expand_range(min_cum, max_cum);
    let mut chart = ChartBuilder::on(&areas[2])
      .caption("Cumulative Error", (CHART_FONT_FAMILY, 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().label_style(LABEL_FONT).draw()?;
    chart.draw_series(LineSeries::new(
      x.iter()
        .zip(sp.cumulative.cumulative_error.iter())
        .map(|(&x, &y)| (x, y)),
      GREEN.stroke_width(2),
    ))?;
    chart.draw_series(once(PathElement::new(
      vec![(x_min, 0.0), (x_max, 0.0)],
      BLACK.stroke_width(1),
    )))?;
  }

  {
    let max_peak = sp
      .regional
      .peak_region_errors
      .iter()
      .copied()
      .filter(|v| v.is_finite())
      .fold(0.0_f64, f64::max);
    let max_tail = sp
      .regional
      .tail_region_errors
      .iter()
      .copied()
      .filter(|v| v.is_finite())
      .fold(0.0_f64, f64::max);
    let max_y = max_peak.max(max_tail);
    let (y_min, y_max) = expand_range(0.0, max_y);
    let mut chart = ChartBuilder::on(&areas[3])
      .caption("Peak & Tail Region Errors", (CHART_FONT_FAMILY, 20))
      .margin(10)
      .x_label_area_size(30)
      .y_label_area_size(50)
      .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
    chart.configure_mesh().label_style(LABEL_FONT).draw()?;
    chart
      .draw_series(LineSeries::new(
        x.iter()
          .zip(sp.regional.peak_region_errors.iter())
          .filter_map(|(&x, &y)| y.is_finite().then_some((x, y))),
        MAGENTA.stroke_width(2),
      ))?
      .label("Peak region")
      .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], MAGENTA.stroke_width(2)));
    chart
      .draw_series(LineSeries::new(
        x.iter()
          .zip(sp.regional.tail_region_errors.iter())
          .filter_map(|(&x, &y)| y.is_finite().then_some((x, y))),
        CYAN.stroke_width(2),
      ))?
      .label("Tail region")
      .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], CYAN.stroke_width(2)));
    chart
      .configure_series_labels()
      .label_font(LABEL_FONT)
      .border_style(BLACK)
      .background_style(WHITE.mix(0.8))
      .draw()?;
  }
  Ok(())
}
