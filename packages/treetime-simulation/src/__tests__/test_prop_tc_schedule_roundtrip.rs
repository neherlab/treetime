#[cfg(test)]
mod tests {
  use generators::{gen_constant, gen_exponential, gen_piecewise};
  use proptest::proptest;
  use treetime_utils::prop_assert_ulps_eq;

  proptest! {
    #[test]
    fn test_prop_tc_schedule_constant_time_after_hazard_roundtrip(schedule in gen_constant(), tau in 0.0..40.0_f64) {
      let actual = schedule.time_after_hazard(0.0, schedule.hazard_between(0.0, tau));
      prop_assert_ulps_eq!(tau, actual, max_ulps = 2);
    }

    #[test]
    fn test_prop_tc_schedule_exponential_time_after_hazard_roundtrip(schedule in gen_exponential(), tau in 0.0..40.0_f64) {
      let actual = schedule.time_after_hazard(0.0, schedule.hazard_between(0.0, tau));
      prop_assert_ulps_eq!(tau, actual, max_ulps = 10);
    }

    #[test]
    fn test_prop_tc_schedule_piecewise_time_after_hazard_roundtrip(schedule in gen_piecewise(), tau in 0.0..40.0_f64) {
      let actual = schedule.time_after_hazard(0.0, schedule.hazard_between(0.0, tau));
      prop_assert_ulps_eq!(tau, actual, max_ulps = 81);
    }
  }

  mod generators {
    use crate::tc_schedule::TcSchedule;
    use itertools::Itertools;
    use proptest::collection::vec;
    use proptest::strategy::Strategy;
    use std::iter::once;

    pub(super) fn gen_constant() -> impl Strategy<Value = TcSchedule> {
      (0.5..20.0_f64).prop_map(|tc| TcSchedule::constant(tc).unwrap())
    }

    pub(super) fn gen_exponential() -> impl Strategy<Value = TcSchedule> {
      (0.5..20.0_f64, 1.0..20.0_f64).prop_map(|(tc0, scale)| TcSchedule::exponential(tc0, scale).unwrap())
    }

    pub(super) fn gen_piecewise() -> impl Strategy<Value = TcSchedule> {
      (0.1..20.0_f64, vec(0.1..10.0_f64, 0..3)).prop_flat_map(|(first, gaps)| {
        let breaks = once(first)
          .chain(gaps)
          .scan(0.0, |end, gap| {
            *end += gap;
            Some(*end)
          })
          .collect_vec();
        let value_count = breaks.len() + 1;
        vec(0.5..20.0_f64, value_count).prop_map(move |values| TcSchedule::piecewise(breaks.clone(), values).unwrap())
      })
    }
  }
}
