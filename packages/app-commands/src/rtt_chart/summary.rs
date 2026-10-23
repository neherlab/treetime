use crate::rtt_chart::data::RttChartData;
use itertools::Itertools;
use treetime::clock::clock_model::{ClockModel, ClockModelStats};
use treetime_utils::cast::f64_to_usize_saturating;
use treetime_utils::datetime::year_fraction::year_fraction_to_date;
use treetime_utils::fmt::float::{float_to_power_of_ten, float_to_power_of_ten_with_error};
use treetime_utils::fmt::string::counted;

const RATE_SIGNIFICANT_DIGITS: usize = 3;

const ROOT_DATE_STEPS_PER_SPAN: f64 = 100.0;

const ROOT_DATE_MAX_DECIMALS: f64 = 3.0;

const ROOT_MONTH_MAX_SPAN: f64 = 10.0;

pub fn root_date(clock_model: &ClockModel) -> Option<f64> {
  let root = clock_model.t_mrca();
  (clock_model.clock_rate() > 0.0 && root.is_finite()).then_some(root)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RttSummary {
  pub rate: Statistic,
  pub root_date: Statistic,
  pub r_squared: Option<Statistic>,
  pub samples: String,
  pub notes: Vec<String>,
}

impl RttSummary {
  pub fn new(data: &RttChartData, clock_model: &ClockModel) -> Self {
    let (first, last) = data.date_range();
    let span = last - first;
    Self {
      rate: rate_statistic(clock_model),
      root_date: root_statistic(root_date(clock_model), span),
      r_squared: clock_model.r_squared().map(|r_squared| Statistic {
        label: "R²".to_owned(),
        value: format!("{r_squared:.3}"),
        detail: None,
      }),
      samples: samples_text(data),
      notes: (data.undated > 0)
        .then(|| {
          let verb = if data.undated == 1 { "is" } else { "are" };
          format!("{} without a date {verb} not shown", counted(data.undated, "sample"))
        })
        .into_iter()
        .collect(),
    }
  }

  pub fn statistics(&self) -> impl Iterator<Item = &Statistic> {
    [&self.rate, &self.root_date].into_iter().chain(&self.r_squared)
  }

  pub fn description(&self) -> String {
    self
      .statistics()
      .map(Statistic::text)
      .chain([self.samples.clone()])
      .chain(self.notes.iter().cloned())
      .join(". ")
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Statistic {
  pub label: String,
  pub value: String,
  pub detail: Option<String>,
}

impl Statistic {
  pub fn text(&self) -> String {
    format!("{} {}", self.label, self.value_text())
  }

  pub fn value_text(&self) -> String {
    match &self.detail {
      Some(detail) => format!("{} ({detail})", self.value),
      None => self.value.clone(),
    }
  }
}

fn rate_statistic(clock_model: &ClockModel) -> Statistic {
  let rate = clock_model.clock_rate();
  let label = match clock_model.stats() {
    ClockModelStats::Estimated(_) => "Clock rate",
    ClockModelStats::Fixed => "Clock rate (fixed)",
  };
  Statistic {
    label: label.to_owned(),
    value: clock_model.rate_std().map_or_else(
      || float_to_power_of_ten(rate, RATE_SIGNIFICANT_DIGITS),
      |rate_std| float_to_power_of_ten_with_error(rate, rate_std, RATE_SIGNIFICANT_DIGITS),
    ),
    detail: Some("substitutions per site per year".to_owned()),
  }
}

fn root_statistic(root: Option<f64>, span: f64) -> Statistic {
  let (value, detail) = match root {
    Some(root) => {
      let decimals = root_date_decimals(span);
      (format!("{root:.decimals$}"), root_month(root, span))
    },
    None => ("None".to_owned(), Some("no positive clock rate".to_owned())),
  };
  Statistic {
    label: "Root date".to_owned(),
    value,
    detail,
  }
}

fn root_date_decimals(span: f64) -> usize {
  if span > 0.0 {
    let decimals = (ROOT_DATE_STEPS_PER_SPAN / span).log10().ceil();
    f64_to_usize_saturating(decimals.clamp(0.0, ROOT_DATE_MAX_DECIMALS))
  } else {
    2
  }
}

fn root_month(root: f64, span: f64) -> Option<String> {
  (root >= 1.0 && span < ROOT_MONTH_MAX_SPAN).then(|| year_fraction_to_date(root).format("%B %Y").to_string())
}

fn samples_text(data: &RttChartData) -> String {
  let fitted = format!("{} in the fit", counted(data.fitted.len(), "sample"));
  if data.excluded.is_empty() {
    fitted
  } else {
    format!("{fitted}, {} excluded by the clock filter", data.excluded.len())
  }
}
