use eyre::Report;
use itertools::{Either, Itertools};
use treetime::clock::rtt::ClockRegressionResult;
use treetime_utils::make_error;

#[derive(Clone, Debug, PartialEq)]
pub struct RttChartData {
  pub fitted: Vec<RttPoint>,
  pub excluded: Vec<RttPoint>,
  pub undated: usize,
}

impl RttChartData {
  pub fn new(results: &[ClockRegressionResult]) -> Result<Self, Report> {
    let tips = results.iter().filter(|result| result.is_leaf).collect_vec();
    let undated = tips.iter().filter(|tip| tip.date.is_none()).count();
    let (excluded, fitted): (Vec<_>, Vec<_>) = tips
      .iter()
      .filter_map(|tip| tip.date.map(|date| (tip.is_outlier, RttPoint { date, div: tip.div })))
      .partition_map(|(is_outlier, point)| {
        if is_outlier {
          Either::Left(point)
        } else {
          Either::Right(point)
        }
      });
    if fitted.is_empty() && excluded.is_empty() {
      return make_error!("the root-to-tip chart needs at least one sample with a date");
    }
    Ok(Self {
      fitted,
      excluded,
      undated,
    })
  }

  pub fn dated(&self) -> impl Iterator<Item = &RttPoint> {
    self.fitted.iter().chain(&self.excluded)
  }

  pub fn date_range(&self) -> (f64, f64) {
    date_extent(self.dated())
  }

  pub fn fit_range(&self) -> (f64, f64) {
    if self.fitted.is_empty() {
      self.date_range()
    } else {
      date_extent(self.fitted.iter())
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RttPoint {
  pub date: f64,
  pub div: f64,
}

fn date_extent<'a>(points: impl Iterator<Item = &'a RttPoint>) -> (f64, f64) {
  points.fold((f64::INFINITY, f64::NEG_INFINITY), |(min, max), point| {
    (min.min(point.date), max.max(point.date))
  })
}
