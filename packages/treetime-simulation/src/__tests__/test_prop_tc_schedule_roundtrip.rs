#[cfg(test)]
mod tests {
  use generators::gen_schedule;
  use proptest::proptest;
  use treetime_utils::prop_assert_ulps_eq;

  proptest! {
    #[test]
    fn test_prop_tc_schedule_inverse_cumulative_hazard_roundtrip(schedule in gen_schedule(), tau in 0.0..40.0_f64) {
      let actual = schedule.inverse_cumulative_hazard(schedule.cumulative_hazard(tau));
      prop_assert_ulps_eq!(tau, actual, max_ulps = 16);
    }
  }

  mod generators {
    use crate::tc_schedule::TcSchedule;
    use itertools::Itertools;
    use proptest::collection::vec;
    use proptest::prop_oneof;
    use proptest::strategy::Strategy;

    pub(super) fn gen_schedule() -> impl Strategy<Value = TcSchedule> {
      prop_oneof![
        (0.5..20.0_f64).prop_map(|tc| TcSchedule::constant(tc).unwrap()),
        (0.5..20.0_f64, 1.0..20.0_f64).prop_map(|(tc0, scale)| TcSchedule::exponential(tc0, scale).unwrap()),
        vec(0.1..20.0_f64, 1..4).prop_flat_map(|breaks| {
          let breaks = breaks.into_iter().sorted_by(f64::total_cmp).dedup().collect_vec();
          let value_count = breaks.len() + 1;
          vec(0.5..20.0_f64, value_count).prop_map(move |values| TcSchedule::piecewise(breaks.clone(), values).unwrap())
        }),
      ]
    }
  }
}
