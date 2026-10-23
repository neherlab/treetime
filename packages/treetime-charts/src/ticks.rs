use eyre::Report;
use num_traits::ToPrimitive;
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating, u32_to_usize};
use treetime_utils::datetime::year_fraction::year_fraction_at_month_start;
use treetime_utils::make_report;

const NICE_MANTISSAS: [(f64, f64); 4] = [(1.0, 0.0), (2.0, 0.0), (2.5, 1.0), (5.0, 0.0)];

const MONTH_STEPS: [u32; 4] = [1, 2, 3, 6];

const YEAR_MANTISSAS: [f64; 3] = [1.0, 2.0, 5.0];

const YEAR_STEP_DECADES: i32 = 8;

const MONTH_NAMES: [&str; 12] = [
  "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
];

const MINUS: &str = "\u{2212}";

pub fn linear_ticks(min: f64, max: f64, max_count: usize) -> Vec<Tick> {
  if !is_range(min, max) || max_count == 0 {
    return vec![];
  }
  let rough_step = (max - min) / max_count.to_f64_exact();
  let base_exponent = rough_step.log10().floor();
  (0..3)
    .flat_map(|power| {
      let exponent = base_exponent + f64::from(power);
      NICE_MANTISSAS.map(|(mantissa, extra_decimals)| (mantissa, extra_decimals, exponent))
    })
    .find_map(|(mantissa, extra_decimals, exponent)| {
      let step = mantissa * 10_f64.powf(exponent);
      let (first, last) = ((min / step).ceil(), (max / step).floor());
      (last - first + 1.0 <= max_count.to_f64_exact()).then(|| {
        let decimals = f64_to_usize_saturating(extra_decimals - exponent);
        multiples(first, last)
          .map(|multiple| Tick {
            value: multiple * step,
            label: decimal_label(multiple * step, multiple == 0.0, decimals),
          })
          .filter(|tick| (min..=max).contains(&tick.value))
          .collect()
      })
    })
    .unwrap_or_default()
}

pub fn date_ticks(min: f64, max: f64, max_count: usize) -> Result<Vec<Tick>, Report> {
  if !is_range(min, max) || max_count == 0 {
    return Ok(vec![]);
  }
  let max_count_f64 = max_count.to_f64_exact();
  for step in MONTH_STEPS {
    if (max - min) * 12.0 / f64::from(step) > max_count_f64 + 1.0 {
      continue;
    }
    let ticks = month_ticks(min, max, step)?;
    if ticks.len() < 2 {
      return Ok(linear_ticks(min, max, max_count));
    }
    if ticks.len() <= max_count {
      return Ok(ticks);
    }
  }
  let year_step = (0..YEAR_STEP_DECADES)
    .flat_map(|decade| YEAR_MANTISSAS.map(|mantissa| mantissa * 10_f64.powi(decade)))
    .find(|&step| (max / step).floor() - (min / step).ceil() + 1.0 <= max_count_f64);
  Ok(year_step.map_or_else(Vec::new, |step| year_ticks(min, max, step)))
}

#[derive(Clone, Debug, PartialEq)]
pub struct Tick {
  pub value: f64,
  pub label: String,
}

fn is_range(min: f64, max: f64) -> bool {
  min.is_finite() && max.is_finite() && min < max
}

fn multiples(first: f64, last: f64) -> impl Iterator<Item = f64> {
  let count = f64_to_usize_saturating(last - first + 1.0);
  (0..count).map(move |index| first + index.to_f64_exact())
}

fn decimal_label(value: f64, is_zero: bool, decimals: usize) -> String {
  if is_zero {
    return "0".to_owned();
  }
  format!("{value:.decimals$}").replace('-', MINUS)
}

fn month_ticks(min: f64, max: f64, step: u32) -> Result<Vec<Tick>, Report> {
  let (first_year, last_year) = (calendar_year(min)?, calendar_year(max)?);
  let months = (first_year..=last_year)
    .flat_map(|year| (1..=12).step_by(u32_to_usize(step)).map(move |month| (year, month)))
    .map(|(year, month)| Ok((year, month, year_fraction_at_month_start(year, month)?)))
    .collect::<Result<Vec<_>, Report>>()?;
  Ok(
    months
      .into_iter()
      .filter(|&(_, _, position)| (min..=max).contains(&position))
      .enumerate()
      .map(|(index, (year, month, position))| Tick {
        value: position,
        label: month_label(year, month, index == 0),
      })
      .collect(),
  )
}

fn calendar_year(value: f64) -> Result<i32, Report> {
  value
    .floor()
    .to_i32()
    .ok_or_else(|| make_report!("Date {value} is outside the calendar range of the chart axis"))
}

fn month_label(year: i32, month: u32, is_first: bool) -> String {
  let name = MONTH_NAMES[u32_to_usize(month - 1)];
  match (month, is_first) {
    (1, _) => year_label(year),
    (_, true) => format!("{name} {}", year_label(year)),
    (_, false) => name.to_owned(),
  }
}

fn year_ticks(min: f64, max: f64, step: f64) -> Vec<Tick> {
  multiples((min / step).ceil(), (max / step).floor())
    .map(|multiple| {
      let value = multiple * step;
      Tick {
        value,
        label: decimal_label(value, multiple == 0.0, 0),
      }
    })
    .filter(|tick| (min..=max).contains(&tick.value))
    .collect()
}

fn year_label(year: i32) -> String {
  year.to_string().replace('-', MINUS)
}
