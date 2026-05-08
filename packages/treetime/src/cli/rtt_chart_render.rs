use crate::clock::clock_model::ClockModel;
use crate::clock::rtt::ClockRegressionResult;
use eyre::Report;
use plotters::coord::Shift;
use plotters::coord::types::RangedCoordf32;
use plotters::prelude::*;
use treetime_utils::fmt::float::float_to_significant_digits;

use crate::cli::rtt_chart::{PointsResult, gather_points};

pub fn draw_chart<'a, DB>(
  results: &[ClockRegressionResult],
  clock_model: &ClockModel,
  drawing_area: &'a DrawingArea<DB, Shift>,
) -> Result<(), Report>
where
  DB: DrawingBackend + 'a,
  <DB as DrawingBackend>::ErrorType: 'static,
{
  let PointsResult {
    norm_points,
    outlier_points,
    line,
    x_min,
    x_max,
    y_min,
    y_max,
  } = gather_points(results, clock_model)?;

  drawing_area.fill(&WHITE)?;
  drawing_area.margin(10, 10, 10, 10);

  let mut chart = ChartBuilder::on(drawing_area)
    .caption(
      "Clock regression",
      ("sans-serif", 40).into_font().style(FontStyle::Bold),
    )
    .x_label_area_size(60)
    .y_label_area_size(60)
    .margin(10)
    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

  chart
    .configure_mesh()
    .light_line_style(WHITE)
    .bold_line_style(RGBAColor(128, 128, 128, 0.5))
    .x_desc("Date")
    .x_label_formatter(&|&date| float_to_significant_digits(date, 6))
    .y_desc("Div")
    .y_label_formatter(&|&div| float_to_significant_digits(div, 3))
    .axis_desc_style(("sans-serif", 20).into_font().style(FontStyle::Bold))
    .x_label_style(("sans-serif", 16))
    .y_label_style(("sans-serif", 16))
    .draw()?;

  let norm_point_color = RGBColor(8, 232, 140);
  chart
    .draw_series(PointSeries::of_element(
      norm_points,
      4,
      &norm_point_color,
      &|c, s, st| {
        EmptyElement::at(c)
          + Circle::new((0, 0), s, st.filled())
          + Circle::new((0, 0), s, RGBAColor(128, 128, 128, 0.25).stroke_width(1))
      },
    ))?
    .label("Samples (norm)")
    .legend(move |(x, y)| {
      Circle::new(
        (x + 10, y),
        4,
        ShapeStyle {
          color: norm_point_color.to_rgba(),
          filled: true,
          stroke_width: 0,
        },
      )
    });

  let outlier_point_color = RGBColor(255, 105, 97);
  chart
    .draw_series(PointSeries::of_element(
      outlier_points,
      4,
      &outlier_point_color,
      &|c, s, st| {
        EmptyElement::at(c)
          + Circle::new((0, 0), s, st.filled())
          + Circle::new((0, 0), s, RGBAColor(128, 128, 128, 0.25).stroke_width(1))
      },
    ))?
    .label("Samples (outliers)")
    .legend(move |(x, y)| {
      Circle::new(
        (x + 10, y),
        4,
        ShapeStyle {
          color: outlier_point_color.to_rgba(),
          filled: true,
          stroke_width: 0,
        },
      )
    });

  let line_color = RGBColor(8, 140, 232);
  chart
    .draw_series(LineSeries::new(
      line,
      ShapeStyle {
        color: line_color.to_rgba(),
        filled: true,
        stroke_width: 3,
      },
    ))?
    .label(format!("Clock regression: {}", clock_model.equation_str()))
    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], line_color));

  add_legend(
    &mut chart,
    format!("tMRCA: {:}", float_to_significant_digits(clock_model.t_mrca(), 5)),
  )?;
  add_legend(
    &mut chart,
    format!("Rate: {:}", float_to_significant_digits(clock_model.clock_rate(), 3)),
  )?;
  add_legend(
    &mut chart,
    format!(
      "Intercept: {:}",
      float_to_significant_digits(clock_model.intercept(), 4)
    ),
  )?;
  if let Some(r_val) = clock_model.r_val() {
    add_legend(&mut chart, format!("R: {:}", float_to_significant_digits(r_val, 4)))?;
    add_legend(
      &mut chart,
      format!("R²: {:}", float_to_significant_digits(r_val.powf(2.0), 4)),
    )?;
  }
  if let Some(chisq) = clock_model.chisq() {
    add_legend(&mut chart, format!("χ²: {chisq:.3e}"))?;
  }

  chart
    .configure_series_labels()
    .label_font(("sans-serif", 20).into_font())
    .margin(5)
    .border_style(RGBAColor(128, 128, 128, 0.5))
    .background_style(RGBAColor(200, 200, 200, 0.5))
    .position(SeriesLabelPosition::UpperLeft)
    .draw()?;

  Ok(())
}

fn add_legend<'a, DB>(
  chart: &mut ChartContext<DB, Cartesian2d<RangedCoordf32, RangedCoordf32>>,
  label: impl Into<String>,
) -> Result<(), Report>
where
  DB: DrawingBackend + 'a,
  <DB as DrawingBackend>::ErrorType: 'static,
{
  chart.draw_series(empty_line())?.label(label);
  Ok(())
}

fn empty_line<DB>() -> LineSeries<DB, (f32, f32)>
where
  DB: DrawingBackend,
{
  LineSeries::<DB, _>::new(Vec::<(f32, f32)>::new(), TRANSPARENT)
}
