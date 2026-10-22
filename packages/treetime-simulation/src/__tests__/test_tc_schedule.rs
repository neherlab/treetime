#[cfg(test)]
mod tests {
  use crate::tc_schedule::TcSchedule;
  use eyre::Report;
  use rstest::rstest;
  use std::f64::consts::{E, LN_10};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::constant_at_zero(      TcSchedule::constant(10.0),                                   0.0, 0.0)]
  #[case::constant(              TcSchedule::constant(10.0),                                   3.0, 0.3)]
  #[case::constant_far(          TcSchedule::constant(10.0),                                  25.0, 2.5)]
  #[case::exponential_one_scale( TcSchedule::exponential(10.0, 5.0),                           5.0, 0.5 * (E - 1.0))]
  #[case::exponential_two_scales(TcSchedule::exponential(10.0, 5.0),                          10.0, 0.5 * (E * E - 1.0))]
  #[case::piecewise_first(       TcSchedule::piecewise(vec![5.0, 6.0], vec![10.0, 1.0, 3.0]),  2.0, 0.2)]
  #[case::piecewise_second(      TcSchedule::piecewise(vec![5.0, 6.0], vec![10.0, 1.0, 3.0]),  5.5, 1.0)]
  #[case::piecewise_last(        TcSchedule::piecewise(vec![5.0, 6.0], vec![10.0, 1.0, 3.0]),  9.0, 2.5)]
  #[trace]
  fn test_tc_schedule_hazard_from_zero_matches_closed_form(
    #[case] schedule: Result<TcSchedule, Report>,
    #[case] tau: f64,
    #[case] expected: f64,
  ) {
    let actual = schedule.unwrap().hazard_between(0.0, tau);
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 1);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::first_segment( 4.0, 10.0)]
  #[case::first_break(   5.0,  1.0)]
  #[case::second_segment(5.5,  1.0)]
  #[case::second_break(  6.0,  3.0)]
  #[case::last_segment(  9.0,  3.0)]
  #[trace]
  fn test_tc_schedule_piecewise_tc_is_right_continuous(#[case] tau: f64, #[case] expected: f64) {
    let schedule = TcSchedule::piecewise(vec![5.0, 6.0], vec![10.0, 1.0, 3.0]).unwrap();
    let actual = schedule.tc(tau);
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 0);
  }

  #[test]
  fn test_tc_schedule_exponential_tc_decays_with_scale() {
    let schedule = TcSchedule::exponential(10.0, 5.0).unwrap();
    let expected = 10.0 / E;
    let actual = schedule.tc(5.0);
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 1);
  }

  #[test]
  fn test_tc_schedule_exponential_hazard_over_no_time_is_zero_when_scale_over_tc0_overflows() {
    let schedule = TcSchedule::exponential(1e-300, 1e300).unwrap();
    let expected = 0.0;
    let actual = schedule.hazard_between(0.0, 0.0);
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 0);
  }

  #[test]
  fn test_tc_schedule_exponential_hazard_is_finite_when_scale_over_tc0_overflows() {
    let schedule = TcSchedule::exponential(1e-300, 1e300).unwrap();
    let expected = 1.0 / 1e-300;
    let actual = schedule.hazard_between(0.0, 1.0);
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 0);
  }

  #[test]
  fn test_tc_schedule_exponential_time_after_hazard_is_finite_when_hazard_times_tc0_overflows() {
    let schedule = TcSchedule::exponential(1e308, 1e308).unwrap();
    let expected = 1e308 * 3.0_f64.ln();
    let actual = schedule.time_after_hazard(0.0, 2.0).unwrap();
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 4);
  }

  #[test]
  fn test_tc_schedule_exponential_hazard_is_finite_when_the_rate_at_the_start_overflows() {
    let schedule = TcSchedule::exponential(1e300, 1.0).unwrap();
    let expected = 750.0 + (E - 1.0).ln() - 300.0 * LN_10;
    let actual = schedule.hazard_between(750.0, 751.0).ln();
    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
  }

  #[test]
  fn test_tc_schedule_exponential_time_after_hazard_is_finite_when_the_rate_at_the_start_overflows() {
    let schedule = TcSchedule::exponential(1e300, 1.0).unwrap();
    let expected = 751.0;
    let actual = schedule
      .time_after_hazard(750.0, schedule.hazard_between(750.0, 751.0))
      .unwrap();
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 1);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::first_break_from_below( 5.0, -1)]
  #[case::first_break_from_above( 5.0,  1)]
  #[case::second_break_from_below(6.0, -1)]
  #[case::second_break_from_above(6.0,  1)]
  #[trace]
  fn test_tc_schedule_hazard_from_zero_is_continuous_at_breaks(#[case] tau: f64, #[case] side: i32) {
    let schedule = TcSchedule::piecewise(vec![5.0, 6.0], vec![10.0, 1.0, 3.0]).unwrap();
    let neighbor = if side < 0 { tau.next_down() } else { tau.next_up() };
    let at_break = schedule.hazard_between(0.0, tau);
    let next_to_break = schedule.hazard_between(0.0, neighbor);
    pretty_assert_abs_diff_eq!(at_break, next_to_break, epsilon = 1e-15);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::below_first_break(      0.0,                           5.0_f64.next_down(),          0)]
  #[case::at_first_break(         0.0,                           5.0,                          1)]
  #[case::at_second_break(        0.0,                           5.0_f64.next_up(),            2)]
  #[case::above_second_break(     0.0,                           5.0_f64.next_up().next_up(),  3)]
  #[case::below_third_break(      0.0,                           (5.0 + 1e-12_f64).next_down(), 12)]
  #[case::at_third_break(         0.0,                           5.0 + 1e-12,                  51)]
  #[case::above_third_break(      0.0,                           (5.0 + 1e-12_f64).next_up(),  76)]
  #[case::from_first_break(       5.0,                           6.0,                          1)]
  #[case::from_second_break(      5.0_f64.next_up(),             6.0,                          1)]
  #[case::from_below_third_break( (5.0 + 1e-12_f64).next_down(), 5.0 + 1e-12,                  0)]
  #[case::across_all_breaks(      5.0_f64.next_down(),           (5.0 + 1e-12_f64).next_up(),  0)]
  #[trace]
  fn test_tc_schedule_piecewise_time_after_hazard_roundtrip_next_to_close_breaks(
    #[case] from: f64,
    #[case] to: f64,
    #[case] max_ulps: u32,
  ) {
    let schedule = TcSchedule::piecewise(vec![5.0, 5.0_f64.next_up(), 5.0 + 1e-12], vec![0.5, 20.0, 0.5, 20.0]).unwrap();
    let actual = schedule.time_after_hazard(from, schedule.hazard_between(from, to)).unwrap();
    pretty_assert_ulps_eq!(to, actual, max_ulps = max_ulps);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::without_redundant_breaks(vec![1.0],           vec![1.0, 2.0_f64.powi(54)])]
  #[case::with_redundant_breaks(   vec![1.0, 2.0, 3.0], vec![1.0, 2.0_f64.powi(54), 2.0_f64.powi(54), 2.0_f64.powi(54)])]
  #[trace]
  fn test_tc_schedule_piecewise_time_after_hazard_is_unchanged_by_breaks_that_keep_tc(
    #[case] breaks: Vec<f64>,
    #[case] values: Vec<f64>,
  ) {
    let schedule = TcSchedule::piecewise(breaks, values).unwrap();
    let expected = 5.0;
    let actual = schedule.time_after_hazard(0.0, 1.0_f64.next_up()).unwrap();
    pretty_assert_ulps_eq!(expected, actual, max_ulps = 0);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::piecewise_negative(  TcSchedule::constant(1.0),         -1.0,     "Hazard after time 2 must be non-negative, got -1")]
  #[case::piecewise_nan(       TcSchedule::constant(1.0),         f64::NAN, "Hazard after time 2 must be non-negative, got NaN")]
  #[case::exponential_negative(TcSchedule::exponential(1.0, 1.0), -1.0,     "Hazard after time 2 must be non-negative, got -1")]
  #[case::exponential_nan(     TcSchedule::exponential(1.0, 1.0), f64::NAN, "Hazard after time 2 must be non-negative, got NaN")]
  #[trace]
  fn test_tc_schedule_time_after_hazard_rejects_invalid_hazard(
    #[case] schedule: Result<TcSchedule, Report>,
    #[case] hazard: f64,
    #[case] expected: &str,
  ) {
    let result = schedule.unwrap().time_after_hazard(2.0, hazard);
    assert_error!(result, expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_tc(              TcSchedule::constant(0.0),                                      "Coalescent time scale Tc must be finite and positive, got 0")]
  #[case::negative_tc(          TcSchedule::constant(-1.0),                                     "Coalescent time scale Tc must be finite and positive, got -1")]
  #[case::infinite_tc(          TcSchedule::constant(f64::INFINITY),                            "Coalescent time scale Tc must be finite and positive, got inf")]
  #[case::nan_piecewise_value(  TcSchedule::piecewise(vec![1.0], vec![1.0, f64::NAN]),         "Coalescent time scale Tc must be finite and positive, got NaN")]
  #[case::zero_exponential_tc0( TcSchedule::exponential(0.0, 1.0),                              "Coalescent time scale Tc must be finite and positive, got 0")]
  #[case::zero_scale(           TcSchedule::exponential(1.0, 0.0),                              "Time scale of an exponential Tc schedule must be finite and positive, got 0")]
  #[case::infinite_scale(       TcSchedule::exponential(1.0, f64::INFINITY),                    "Time scale of an exponential Tc schedule must be finite and positive, got inf")]
  #[case::unsorted_breaks(      TcSchedule::piecewise(vec![6.0, 5.0], vec![1.0, 2.0, 3.0]),     "Breaks of a piecewise Tc schedule must be strictly increasing, got [6.0, 5.0]")]
  #[case::repeated_break(       TcSchedule::piecewise(vec![5.0, 5.0], vec![1.0, 2.0, 3.0]),     "Breaks of a piecewise Tc schedule must be strictly increasing, got [5.0, 5.0]")]
  #[case::nan_break(            TcSchedule::piecewise(vec![f64::NAN], vec![1.0, 2.0]),          "Breaks of a piecewise Tc schedule must be finite, got NaN")]
  #[case::infinite_break(       TcSchedule::piecewise(vec![f64::INFINITY], vec![1.0, 2.0]),     "Breaks of a piecewise Tc schedule must be finite, got inf")]
  #[case::too_few_values(       TcSchedule::piecewise(vec![5.0], vec![10.0]),                   "A piecewise Tc schedule with 1 break needs 2 values, got 1")]
  #[case::too_many_values(      TcSchedule::piecewise(vec![5.0], vec![10.0, 1.0, 3.0]),         "A piecewise Tc schedule with 1 break needs 2 values, got 3")]
  #[case::no_values(            TcSchedule::piecewise(vec![], vec![]),                          "A piecewise Tc schedule with 0 breaks needs 1 value, got 0")]
  #[trace]
  fn test_tc_schedule_rejects_invalid_input(#[case] schedule: Result<TcSchedule, Report>, #[case] expected: &str) {
    assert_error!(schedule, expected);
  }
}
