use itertools::Itertools;
use ndarray::Array1;

pub fn combined_range(arrays: &[&Array1<f64>]) -> (f64, f64) {
  arrays
    .iter()
    .flat_map(|array| array.iter())
    .minmax()
    .into_option()
    .map_or((0.0, 1.0), |(min, max)| (*min, *max))
}

pub fn expand_range(min_value: f64, max_value: f64) -> (f64, f64) {
  if !min_value.is_finite() || !max_value.is_finite() {
    return (-1.0, 1.0);
  }
  if max_value <= min_value {
    let base = if min_value.abs() > 1.0 { min_value.abs() } else { 1.0 };
    return (min_value - base, max_value + base);
  }
  let span = max_value - min_value;
  let padding = if span * 0.05 < 1e-9 { 1e-9 } else { span * 0.05 };
  (min_value - padding, max_value + padding)
}

pub fn tolerance_label(value: f64) -> String {
  let idx = (value + 0.5).floor() as i32;
  match idx {
    0 => "Strict".to_owned(),
    1 => "Moderate".to_owned(),
    2 => "Loose".to_owned(),
    _ => format!("{value:.1}"),
  }
}
