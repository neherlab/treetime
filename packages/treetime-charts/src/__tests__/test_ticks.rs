#[cfg(test)]
mod tests {
  use crate::ticks::{Tick, date_ticks, linear_ticks};
  use helpers::labels;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::small_divergences((0.0,   0.00047, 6), vec!["0", "0.0001", "0.0002", "0.0003", "0.0004"])]
  #[case::offset_range(     (0.11,  0.175,   6), vec!["0.12", "0.14", "0.16"])]
  #[case::quarter_steps(    (0.0,   1.0,     5), vec!["0", "0.25", "0.50", "0.75", "1.00"])]
  #[case::crosses_zero(     (-0.02, 0.05,    8), vec!["\u{2212}0.02", "\u{2212}0.01", "0", "0.01", "0.02", "0.03", "0.04", "0.05"])]
  #[case::integers(         (0.0,   120.0,   4), vec!["0", "50", "100"])]
  #[trace]
  fn test_linear_ticks_use_one_number_format_per_axis(
    #[case] (min, max, max_count): (f64, f64, usize),
    #[case] expected: Vec<&str>,
  ) {
    assert_eq!(expected, labels(&linear_ticks(min, max, max_count)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::empty_range( (1.0,      1.0, 5))]
  #[case::inverted(    (2.0,      1.0, 5))]
  #[case::not_finite(  (f64::NAN, 1.0, 5))]
  #[case::no_ticks(    (0.0,      1.0, 0))]
  #[trace]
  fn test_linear_ticks_of_a_degenerate_axis_are_empty(#[case] (min, max, max_count): (f64, f64, usize)) {
    assert_eq!(Vec::<Tick>::new(), linear_ticks(min, max, max_count));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::quarters_with_the_year_at_january(
    (2014.13, 2016.33, 10),
    vec!["Apr 2014", "Jul", "Oct", "2015", "Apr", "Jul", "Oct", "2016", "Apr"],
  )]
  #[case::first_tick_names_its_year(
    (2015.3, 2015.9, 10),
    vec!["May 2015", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"],
  )]
  #[case::half_years(
    (2013.7, 2016.6, 8),
    vec!["2014", "Jul", "2015", "Jul", "2016", "Jul"],
  )]
  #[case::decades(
    (1940.0, 2020.5, 10),
    vec!["1940", "1950", "1960", "1970", "1980", "1990", "2000", "2010", "2020"],
  )]
  #[case::years_around_the_era(
    (-12.0, 7.0, 5),
    vec!["\u{2212}10", "\u{2212}5", "0", "5"],
  )]
  #[case::shorter_than_two_months_uses_decimal_years(
    (2019.999, 2020.042, 10),
    vec!["2020.000", "2020.005", "2020.010", "2020.015", "2020.020", "2020.025", "2020.030", "2020.035", "2020.040"],
  )]
  #[trace]
  fn test_date_ticks_label_calendar_boundaries(
    #[case] (min, max, max_count): (f64, f64, usize),
    #[case] expected: Vec<&str>,
  ) {
    assert_eq!(expected, labels(&date_ticks(min, max, max_count).unwrap()));
  }

  #[test]
  fn test_date_ticks_sit_at_the_first_instant_of_each_month() {
    let positions = date_ticks(2014.13, 2014.8, 4)
      .unwrap()
      .into_iter()
      .map(|tick| tick.value)
      .collect_vec();
    let expected = vec![
      2014.0 + f64::from(31 + 28) / 365.0,
      2014.0 + f64::from(31 + 28 + 31 + 30) / 365.0,
      2014.0 + f64::from(31 + 28 + 31 + 30 + 31 + 30) / 365.0,
      2014.0 + f64::from(31 + 28 + 31 + 30 + 31 + 30 + 31 + 31) / 365.0,
    ];
    assert_eq!(expected, positions);
  }

  proptest! {
    #[test]
    fn test_prop_linear_ticks_lie_in_the_range_in_order_with_one_format(
      min in -1e3_f64..1e3,
      span_exponent in -8_f64..4.0,
      max_count in 3_usize..12,
    ) {
      let max = min + 10_f64.powf(span_exponent);
      prop_assume!(max > min);
      let ticks = linear_ticks(min, max, max_count);

      prop_assert!((1..=max_count).contains(&ticks.len()), "{ticks:?}");
      prop_assert!(ticks.iter().all(|tick| (min..=max).contains(&tick.value)), "{ticks:?}");
      prop_assert!(ticks.iter().tuple_windows().all(|(a, b)| a.value < b.value), "{ticks:?}");
      let decimals = ticks
        .iter()
        .filter(|tick| tick.label != "0")
        .map(|tick| tick.label.split_once('.').map_or(0, |(_, fraction)| fraction.len()))
        .unique()
        .collect_vec();
      prop_assert!(decimals.len() <= 1, "{ticks:?}");
    }

    #[test]
    fn test_prop_date_ticks_lie_in_the_range_in_order(
      min in 1000_f64..3000.0,
      span_exponent in -1.5_f64..2.5,
      max_count in 3_usize..12,
    ) {
      let max = min + 10_f64.powf(span_exponent);
      let ticks = date_ticks(min, max, max_count).unwrap();

      prop_assert!((1..=max_count).contains(&ticks.len()), "{ticks:?}");
      prop_assert!(ticks.iter().all(|tick| (min..=max).contains(&tick.value)), "{ticks:?}");
      prop_assert!(ticks.iter().tuple_windows().all(|(a, b)| a.value < b.value), "{ticks:?}");
    }
  }

  mod helpers {
    use crate::ticks::Tick;

    pub(super) fn labels(ticks: &[Tick]) -> Vec<&str> {
      ticks.iter().map(|tick| tick.label.as_str()).collect()
    }
  }
}
