use crate::clock::graph_regression::BaseRegressionResult;
use crate::clock::run_clock_model::{ClockModel, RootToTipResult};
use itertools::Itertools;
use log::{info, warn};
use num_traits::clamp;
use rgb::RGB8;
use textplots::{Chart, ColorPlot, Shape};

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
  info!(
    "  distance(t) = {:.4}t {:} {:.4}",
    slope,
    if intercept < &0.0 { "-" } else { "+" },
    intercept.abs()
  );
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
