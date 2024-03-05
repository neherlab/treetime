use crate::commands::clock::graph_regression::BaseRegressionResult;
use crate::commands::clock::run_clock_model::{ClockModel, RootToTipResult};
use eyre::Report;
use itertools::Itertools;
use log::{info, warn};
use num_traits::clamp;
use plotters::{
  backend::SVGBackend,
  chart::{ChartBuilder, SeriesLabelPosition},
  drawing::IntoDrawingArea,
  element::{Circle, EmptyElement, PathElement},
  series::{LineSeries, PointSeries},
  style::{Color, IntoFont, RGBColor, ShapeStyle, BLACK, WHITE},
};
use rgb::RGB8;
use std::path::Path;
use textplots::{Chart, ColorPlot, Shape};

pub fn write_rtt_svg_chart(
  filename: impl AsRef<Path>,
  rtt: &[RootToTipResult],
  clock_model: &ClockModel,
) -> Result<(), Report> {
  let ClockModel {
    regression: BaseRegressionResult {
      slope,
      intercept,
      chisq,
      hessian,
      cov,
    },
    r_val,
  } = clock_model;

  let rtt_points = rtt
    .iter()
    .map(|result| (result.date as f32, result.clock_deviation as f32))
    .collect_vec();

  let (x_min, x_max) = match rtt_points.iter().minmax_by_key(|x| x.0).into_option() {
    None => {
      warn!("When drawing root-to-tip chart: unable to find minimum and maximum of x axis");
      return Ok(());
    }
    Some(((x_min, _), (x_max, _))) => (x_min * 0.999, x_max * 1.001),
  };

  let (y_min, y_max) = match rtt_points.iter().minmax_by_key(|x| x.1).into_option() {
    None => {
      warn!("When drawing root-to-tip chart: unable to find minimum and maximum of y axis");
      return Ok(());
    }
    Some(((_, y_min), (_, y_max))) => (y_min * 0.999, y_max * 1.001),
  };

  let f = |t| *slope as f32 * t + *intercept as f32;
  let rtt_line = [(x_min, f(x_min)), (x_max, f(x_max))];

  let root = {
    let root = SVGBackend::new(filename.as_ref(), (640, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    root.margin(10, 10, 10, 10)
  };

  let mut chart = ChartBuilder::on(&root)
    .caption("Root-to-tip regression", ("sans-serif", 20).into_font())
    .x_label_area_size(30)
    .y_label_area_size(60)
    .margin(5)
    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;

  chart
    .configure_mesh()
    .x_desc("Date")
    .y_desc("Root-to-tip distance")
    .axis_desc_style(("sans-serif", 12))
    .draw()?;

  let line_color = RGBColor(8, 140, 232);
  chart
    .draw_series(LineSeries::new(rtt_line, &line_color))?
    .label(format!(
      "Root-to-tip regression: {}",
      rtt_equation_str(*slope, *intercept)
    ))
    .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], line_color));

  let point_color = RGBColor(255, 105, 97);
  chart
    .draw_series(PointSeries::of_element(rtt_points, 2, &point_color, &|c, s, st| {
      EmptyElement::at(c) + Circle::new((0, 0), s, st.filled())
    }))?
    .label("Actual points")
    .legend(move |(x, y)| {
      Circle::new(
        (x + 10, y),
        2,
        ShapeStyle {
          color: point_color.to_rgba(),
          filled: true,
          stroke_width: 0,
        },
      )
    });

  chart
    .configure_series_labels()
    .border_style(BLACK)
    .background_style(WHITE)
    .position(SeriesLabelPosition::UpperRight)
    .draw()?;

  root.present()?;

  Ok(())
}

pub fn draw_rtt_console_chart(rtt: &[RootToTipResult], clock_model: &ClockModel) {
  let ClockModel {
    regression: BaseRegressionResult {
      slope,
      intercept,
      chisq,
      hessian,
      cov,
    },
    r_val,
  } = clock_model;

  let rtt_points = rtt
    .iter()
    .map(|result| (result.date as f32, result.clock_deviation as f32))
    .collect_vec();

  let (x_min, x_max) = match rtt_points.iter().minmax_by_key(|x| x.0).into_option() {
    None => {
      warn!("When drawing root-to-tip chart: unable to find minimum and maximum of time axis");
      return;
    }
    Some(((x_min, _), (x_max, _))) => (x_min * 0.999, x_max * 1.001),
  };

  let rtt_line = Box::new(|x| *slope as f32 * x + *intercept as f32);

  info!("Root to tip regression:");
  info!("  {}", rtt_equation_str(*slope, *intercept));
  info!("Root date: {:.1}", -intercept / slope);
  info!("Rate:      {:.4}", slope);
  info!("R²:        {:.4}", r_val.powf(2.0));

  let (width, height) = match crossterm::terminal::size() {
    Ok((width, height)) => {
      let width = clamp(width, 32, 1024) as u32;
      let height = clamp(height, 32, 1024) as u32;
      (width, height)
    }
    Err(err) => {
      warn!("When drawing root-to-tip chart: unable to find terminal dimensions: {err}");
      return;
    }
  };

  Chart::new(width, height, x_min, x_max)
    .linecolorplot(&Shape::Points(&rtt_points), RGB8 { r: 255, g: 105, b: 97 })
    .linecolorplot(&Shape::Continuous(rtt_line), RGB8 { r: 8, g: 140, b: 232 })
    .display();
}

fn rtt_equation_str(slope: f64, intercept: f64) -> String {
  format!(
    "distance(t) = {:.4}t {:} {:.4}",
    slope,
    if intercept < 0.0 { "-" } else { "+" },
    intercept.abs()
  )
}
