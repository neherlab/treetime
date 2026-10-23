use crate::rtt_chart::data::RttChartData;
use crate::rtt_chart::summary::root_date;
use treetime::clock::clock_model::ClockModel;

const DATE_PADDING: f64 = 0.04;

const DIVERGENCE_PADDING: f64 = 0.06;

const EMPTY_DATE_SPAN_PADDING: f64 = 0.5;

#[derive(Clone, Debug, PartialEq)]
pub struct RttFrame {
  pub dates: (f64, f64),
  pub divergences: (f64, f64),
  pub fit: [(f64, f64); 2],
  pub extrapolation: Option<[(f64, f64); 2]>,
  pub root: RootPlacement,
}

impl RttFrame {
  pub fn new(data: &RttChartData, clock_model: &ClockModel) -> Self {
    let (first, last) = data.date_range();
    let root = root_placement(root_date(clock_model), first, last);
    let (low, high) = match root {
      RootPlacement::OnChart(root) => (first.min(root), last.max(root)),
      RootPlacement::BeforeChart(_) | RootPlacement::AfterChart(_) | RootPlacement::Undefined => (first, last),
    };
    let dates = padded(low, high, DATE_PADDING, EMPTY_DATE_SPAN_PADDING);

    let line = |date: f64| (date, clock_model.div(date));
    let (fit_first, fit_last) = data.fit_range();
    let fit = [line(fit_first), line(fit_last)];
    let extrapolation = match root {
      RootPlacement::OnChart(root) if root < fit_first => Some([(root, 0.0), line(fit_first)]),
      RootPlacement::OnChart(root) if root > fit_last => Some([line(fit_last), (root, 0.0)]),
      RootPlacement::BeforeChart(_) => Some([line(dates.0), line(fit_first)]),
      RootPlacement::AfterChart(_) => Some([line(fit_last), line(dates.1)]),
      RootPlacement::OnChart(_) | RootPlacement::Undefined => None,
    };

    let line_points = fit.iter().chain(extrapolation.iter().flatten());
    let (div_low, div_high) = data
      .dated()
      .map(|point| point.div)
      .chain(line_points.map(|&(_, div)| div))
      .fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), div| {
        (min.min(div), max.max(div))
      });
    let divergence_pad_floor = div_high.abs().max(div_low.abs()) * DIVERGENCE_PADDING;
    let divergences = padded(div_low, div_high, DIVERGENCE_PADDING, divergence_pad_floor);

    Self {
      dates,
      divergences,
      fit,
      extrapolation,
      root,
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RootPlacement {
  OnChart(f64),
  BeforeChart(f64),
  AfterChart(f64),
  Undefined,
}

fn root_placement(root: Option<f64>, first: f64, last: f64) -> RootPlacement {
  let reach = if last > first { last - first } else { 1.0 };
  match root {
    Some(root) if root < first - reach => RootPlacement::BeforeChart(root),
    Some(root) if root > last + reach => RootPlacement::AfterChart(root),
    Some(root) => RootPlacement::OnChart(root),
    None => RootPlacement::Undefined,
  }
}

fn padded(low: f64, high: f64, fraction: f64, empty_span_padding: f64) -> (f64, f64) {
  let span = high - low;
  let padding = if span > 0.0 {
    span * fraction
  } else {
    empty_span_padding
  };
  let padding = if padding > 0.0 { padding } else { 1.0 };
  (low - padding, high + padding)
}
