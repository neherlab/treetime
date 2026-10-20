#[cfg(test)]
mod tests {
  use crate::coalescent::knot_times::{KnotCursor, KnotGrid, KnotTimes};
  use itertools::Itertools;
  use ordered_float::OrderedFloat;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use std::cell::Cell;
  use treetime_utils::cast::ToF64Exact;

  #[test]
  fn test_knot_times_keep_listed_times_between_grid_runs() {
    let grid = KnotGrid {
      origin: -0.3,
      step: 0.1,
    };
    let on_grid = |index: usize| grid.origin + grid.step * index.to_f64_exact();
    let times = vec![
      -0.45,
      -0.4,
      on_grid(1),
      on_grid(2),
      on_grid(3),
      on_grid(3),
      0.05,
      on_grid(4),
      on_grid(5),
      0.31,
    ];

    let knots = KnotTimes::new(&times, Some(grid));

    assert_eq!(
      helpers::bits(&times),
      helpers::bits(&(0..knots.len()).map(|i| knots.at(i)).collect_vec())
    );
    assert_eq!(helpers::bits(&times), helpers::bits(&knots.iter().collect_vec()));
    assert_eq!((times[0], times[9]), (knots.first(), knots.last()));
  }

  #[test]
  fn test_knot_times_partition_point_over_repeated_knots_evaluates_logarithmically() {
    let count = 2000;
    let knots = helpers::repeated_knots(count);
    let limit = 4 * (usize::BITS - knots.len().leading_zeros());

    let costs = (0..count)
      .map(|k| {
        let x = k.to_f64_exact() + 0.5;
        let calls = Cell::new(0_u32);
        let position = knots.partition_point(|t| {
          calls.set(calls.get() + 1);
          t <= x
        });
        (position, calls.get() <= limit)
      })
      .collect_vec();

    assert_eq!((0..count).map(|k| (2 * k + 2, true)).collect_vec(), costs);
  }

  #[test]
  fn test_knot_times_sorted_sweep_over_repeated_knots_evaluates_linearly() {
    let count = 2000;
    let knots = helpers::repeated_knots(count);
    let calls = Cell::new(0_usize);

    let positions = (0..count)
      .scan(KnotCursor::default(), |cursor, k| {
        let x = k.to_f64_exact() + 0.5;
        *cursor = knots.advance_while(*cursor, |t| {
          calls.set(calls.get() + 1);
          t <= x
        });
        Some(cursor.position())
      })
      .collect_vec();

    assert_eq!((0..count).map(|k| 2 * k + 1).collect_vec(), positions);
    assert!(calls.get() <= 8 * count, "{} evaluations", calls.get());
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn test_prop_knot_times_reproduce_times_and_partition_points(
      origin in -3000.0_f64..3000.0,
      step in 1e-6_f64..10.0,
      indices in prop::collection::vec(0_usize..2000, 1..60),
      offsets in prop::collection::vec(0.0_f64..1.0, 0..20),
      probes in prop::collection::vec(-0.1_f64..1.1, 1..20),
    ) {
      let on_grid = |index: usize| origin + step * index.to_f64_exact();
      let times = indices
        .iter()
        .map(|&index| on_grid(index))
        .chain(offsets.iter().map(|&offset| on_grid(0) + offset * step * 2000.0))
        .sorted_by_key(|&t| OrderedFloat(t))
        .collect_vec();
      let span = times[times.len() - 1] - times[0];

      for grid in [Some(KnotGrid { origin, step }), None] {
        let knots = KnotTimes::new(&times, grid);
        prop_assert_eq!(times.len(), knots.len());
        prop_assert_eq!(helpers::bits(&times), helpers::bits(&knots.iter().collect_vec()));
        for (position, (&a, &b)) in times.iter().tuple_windows().enumerate() {
          let (start, end) = knots.segment(position);
          prop_assert_eq!((a.to_bits(), b.to_bits()), (start.to_bits(), end.to_bits()));
        }
        for &probe in &probes {
          let x = times[0] + probe * span;
          prop_assert_eq!(times.partition_point(|&t| t <= x), knots.partition_point(|t| t <= x));
          prop_assert_eq!(times.partition_point(|&t| t < x), knots.partition_point(|t| t < x));
        }
        let probes = probes.iter().map(|&probe| times[0] + probe * span).sorted_by_key(|&x| OrderedFloat(x));
        let mut cursor = KnotCursor::default();
        for x in probes.filter(|&x| x >= times[0]) {
          cursor = knots.advance_while(cursor, |t| t <= x);
          let position = times.partition_point(|&t| t <= x) - 1;
          prop_assert_eq!(position, cursor.position());
          if position + 1 < times.len() {
            prop_assert_eq!(knots.segment(position), knots.segment_at(cursor));
          }
        }
        for &t in &times {
          prop_assert_eq!(times.partition_point(|&s| s <= t), knots.partition_point(|s| s <= t));
        }
      }
    }
  }

  mod helpers {
    use crate::coalescent::knot_times::{KnotGrid, KnotTimes};
    use itertools::Itertools;
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn repeated_knots(count: usize) -> KnotTimes {
      let grid = KnotGrid { origin: 0.0, step: 1.0 };
      let times = (0..=count).flat_map(|k| [k.to_f64_exact(); 2]).collect_vec();
      KnotTimes::new(&times, Some(grid))
    }

    pub(super) fn bits(times: &[f64]) -> Vec<u64> {
      times.iter().map(|t| t.to_bits()).collect()
    }
  }
}
