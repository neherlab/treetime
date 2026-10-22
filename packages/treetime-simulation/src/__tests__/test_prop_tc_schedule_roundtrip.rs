#[cfg(test)]
mod tests {
  use generators::{gen_constant, gen_exponential, gen_piecewise, gen_span};
  use proptest::proptest;
  use treetime_utils::prop_assert_ulps_eq;

  proptest! {
    #[test]
    fn test_prop_tc_schedule_constant_time_after_hazard_roundtrip(schedule in gen_constant(), span in gen_span()) {
      let actual = schedule.time_after_hazard(span.from, schedule.hazard_between(span.from, span.to)).unwrap();
      prop_assert_ulps_eq!(span.to, actual, max_ulps = 4);
    }

    #[test]
    fn test_prop_tc_schedule_exponential_time_after_hazard_roundtrip(schedule in gen_exponential(), span in gen_span()) {
      let actual = schedule.time_after_hazard(span.from, schedule.hazard_between(span.from, span.to)).unwrap();
      prop_assert_ulps_eq!(span.to, actual, max_ulps = 19);
    }

    #[test]
    fn test_prop_tc_schedule_piecewise_time_after_hazard_roundtrip(schedule in gen_piecewise(), span in gen_span()) {
      let actual = schedule.time_after_hazard(span.from, schedule.hazard_between(span.from, span.to)).unwrap();
      prop_assert_ulps_eq!(span.to, actual, max_ulps = 124);
    }
  }

  mod generators {
    use crate::tc_schedule::TcSchedule;
    use itertools::Itertools;
    use proptest::collection::vec;
    use proptest::prop_oneof;
    use proptest::strategy::Strategy;
    use std::iter::{once, successors};

    pub(super) fn gen_constant() -> impl Strategy<Value = TcSchedule> {
      (0.5..20.0_f64).prop_map(|tc| TcSchedule::constant(tc).unwrap())
    }

    pub(super) fn gen_exponential() -> impl Strategy<Value = TcSchedule> {
      (0.5..20.0_f64, 1.0..20.0_f64).prop_map(|(tc0, scale)| TcSchedule::exponential(tc0, scale).unwrap())
    }

    pub(super) fn gen_piecewise() -> impl Strategy<Value = TcSchedule> {
      (0.1..20.0_f64, vec(gen_gap(), 0..3)).prop_flat_map(|(first, gaps)| {
        let breaks = once(first)
          .chain(gaps.into_iter().scan(first, |last, gap| {
            *last = gap.after(*last);
            Some(*last)
          }))
          .collect_vec();
        let value_count = breaks.len() + 1;
        vec(0.5..20.0_f64, value_count).prop_map(move |values| TcSchedule::piecewise(breaks.clone(), values).unwrap())
      })
    }

    pub(super) fn gen_span() -> impl Strategy<Value = Span> {
      (0.0..40.0_f64, 0.0..40.0_f64).prop_map(|(first, second)| Span {
        from: first.min(second),
        to: first.max(second),
      })
    }

    fn gen_gap() -> impl Strategy<Value = Gap> {
      prop_oneof![
        (0.1..10.0_f64).prop_map(Gap::Width),
        (1e-14..1e-9_f64).prop_map(Gap::Width),
        (1..4_usize).prop_map(Gap::Ulps),
      ]
    }

    #[derive(Clone, Debug)]
    pub(super) struct Span {
      pub(super) from: f64,
      pub(super) to: f64,
    }

    #[derive(Clone, Debug)]
    enum Gap {
      Width(f64),
      Ulps(usize),
    }

    impl Gap {
      fn after(&self, time: f64) -> f64 {
        match *self {
          Self::Width(width) => time + width,
          Self::Ulps(count) => successors(Some(time), |&time| Some(time.next_up())).nth(count).unwrap(),
        }
      }
    }
  }
}
