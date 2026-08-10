#[cfg(test)]
mod tests {
  use crate::timetree::optimization::polytomy::sweep::{Lineage, SubtreePlan, simulate_subtree};
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_utils::sync::random::get_random_number_generator;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  /// A merger rate that ignores time, so tests can reason about the sweep alone.
  fn const_merger_rate(rate: f64) -> PiecewiseConstantFn {
    PiecewiseConstantFn::new(array![], array![rate])
  }

  fn lineage(time: f64, mutations: u32) -> Lineage {
    Lineage { time, mutations }
  }

  /// Every original child must appear exactly once in the forest the plan describes, either
  /// as a root or beneath exactly one merger.
  fn assert_every_child_placed_once(plan: &SubtreePlan, n_children: usize) {
    let mut seen: Vec<usize> = plan.roots.clone();
    for merger in &plan.mergers {
      seen.push(merger.left);
      seen.push(merger.right);
    }

    let unique: BTreeSet<usize> = seen.iter().copied().collect();
    assert_eq!(
      unique.len(),
      seen.len(),
      "a lineage was placed more than once: {seen:?}"
    );

    let total = n_children + plan.mergers.len();
    let expected: BTreeSet<usize> = (0..total).collect();
    assert_eq!(unique, expected, "plan does not account for every lineage exactly once");
  }

  #[test]
  fn test_sweep_two_children_is_noop() -> Result<(), Report> {
    let mut rng = get_random_number_generator(Some(1));
    let children = [lineage(10.0, 0), lineage(9.0, 0)];

    let plan = simulate_subtree(&children, 0.0, 1.0, &const_merger_rate(1.0), &mut rng)?;

    assert!(plan.mergers.is_empty(), "a bifurcation is not a polytomy");
    assert_eq!(plan.roots, vec![0, 1]);
    Ok(())
  }

  #[test]
  fn test_sweep_without_time_window_is_noop() -> Result<(), Report> {
    let mut rng = get_random_number_generator(Some(1));
    // Every child is at or older than the parent, so there is nowhere to put a merger.
    let children = [lineage(5.0, 0), lineage(4.0, 0), lineage(3.0, 0)];

    let plan = simulate_subtree(&children, 5.0, 1.0, &const_merger_rate(1e6), &mut rng)?;

    assert!(plan.mergers.is_empty(), "no window means no resolution");
    assert_eq!(plan.roots, vec![0, 1, 2]);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::parent_nan(       f64::NAN,       1.0, vec![lineage(3.0, 0),       lineage(2.0, 0), lineage(1.0, 0)], "Polytomy parent time must be finite, got NaN")]
  #[case::mutation_negative(0.0,            -1.0, vec![lineage(3.0, 0),       lineage(2.0, 0), lineage(1.0, 0)], "Polytomy mutation rate must be finite and non-negative, got -1")]
  #[case::mutation_nan(     0.0,        f64::NAN, vec![lineage(3.0, 0),       lineage(2.0, 0), lineage(1.0, 0)], "Polytomy mutation rate must be finite and non-negative, got NaN")]
  #[case::child_nan(        0.0,             1.0, vec![lineage(3.0, 0), lineage(f64::NAN, 0), lineage(1.0, 0)], "Polytomy child 1 time must be finite, got NaN")]
  #[trace]
  fn test_sweep_rejects_invalid_inputs(
    #[case] t_stop: f64,
    #[case] mutation_rate: f64,
    #[case] children: Vec<Lineage>,
    #[case] expected: &str,
  ) {
    let mut rng = get_random_number_generator(Some(1));

    assert_error!(
      simulate_subtree(&children, t_stop, mutation_rate, &const_merger_rate(1.0), &mut rng),
      expected
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(-1.0,          "Polytomy merger rate must be finite and non-negative at calendar time 2.500000e0, got -1.000000e0")]
  #[case::nan(     f64::NAN,      "Polytomy merger rate must be finite and non-negative at calendar time 2.500000e0, got NaN")]
  #[case::infinite(f64::INFINITY, "Polytomy merger rate must be finite and non-negative at calendar time 2.500000e0, got inf")]
  #[trace]
  fn test_sweep_rejects_invalid_merger_rates(#[case] rate: f64, #[case] expected: &str) {
    let children = vec![lineage(3.0, 0), lineage(2.0, 0), lineage(1.0, 0)];
    let mut rng = get_random_number_generator(Some(1));

    assert_error!(
      simulate_subtree(&children, 0.0, 1.0, &const_merger_rate(rate), &mut rng),
      expected
    );
  }

  #[test]
  fn test_sweep_rejects_overflowing_component_rate() {
    let children = vec![lineage(3.0, u32::MAX), lineage(3.0, 0), lineage(3.0, 0)];
    let mut rng = get_random_number_generator(Some(1));

    assert_error!(
      simulate_subtree(&children, 0.0, f64::MAX, &const_merger_rate(1.0), &mut rng,),
      "Polytomy event rates must be finite at calendar time 1.500000e0, got mutation rate inf and merger rate 1.000000e0"
    );
  }

  #[test]
  fn test_sweep_rejects_overflowing_total_rate() {
    let children = vec![lineage(3.0, 1), lineage(3.0, 0), lineage(3.0, 0)];
    let component_rate = f64::MAX * 0.75;
    let mut rng = get_random_number_generator(Some(1));

    assert_error!(
      simulate_subtree(
        &children,
        0.0,
        component_rate,
        &const_merger_rate(component_rate),
        &mut rng,
      ),
      "Polytomy total event rate must be finite at calendar time 1.500000e0, got inf"
    );
  }

  #[test]
  fn test_sweep_is_reproducible_under_the_same_seed() -> Result<(), Report> {
    let children: Vec<Lineage> = (0..8_u32).map(|i| lineage(10.0 - f64::from(i), i % 3)).collect();

    let run = |seed: u64| -> Result<SubtreePlan, Report> {
      let mut rng = get_random_number_generator(Some(seed));
      simulate_subtree(&children, -20.0, 0.5, &const_merger_rate(0.4), &mut rng)
    };

    assert_eq!(run(7)?, run(7)?, "same seed must produce the same plan");
    assert_ne!(
      run(7)?,
      run(8)?,
      "different seeds must produce different plans for this configuration"
    );
    Ok(())
  }

  #[test]
  fn test_sweep_preserves_elapsed_events_under_calendar_translation() -> Result<(), Report> {
    let children = vec![
      lineage(10.0, 0),
      lineage(9.5, 2),
      lineage(8.0, 1),
      lineage(7.0, 0),
      lineage(6.0, 3),
      lineage(2.0, 0),
    ];
    let merger_rate = PiecewiseConstantFn::new(array![3.0, 8.0], array![0.2, 0.5, 0.1]);
    let offset = 2048.0;
    let shifted_children: Vec<Lineage> = children
      .iter()
      .map(|child| lineage(child.time + offset, child.mutations))
      .collect();
    let shifted_merger_rate = PiecewiseConstantFn::new(array![3.0 + offset, 8.0 + offset], array![0.2, 0.5, 0.1]);

    let mut base_rng = get_random_number_generator(Some(17));
    let base = simulate_subtree(&children, -5.0, 0.7, &merger_rate, &mut base_rng)?;
    let mut shifted_rng = get_random_number_generator(Some(17));
    let shifted = simulate_subtree(
      &shifted_children,
      -5.0 + offset,
      0.7,
      &shifted_merger_rate,
      &mut shifted_rng,
    )?;

    // Translation oracle: all hazards and branch lengths depend only on time differences.
    let base_pairs: Vec<(usize, usize)> = base.mergers.iter().map(|merger| (merger.left, merger.right)).collect();
    let shifted_pairs: Vec<(usize, usize)> = shifted
      .mergers
      .iter()
      .map(|merger| (merger.left, merger.right))
      .collect();
    assert_eq!(base_pairs, shifted_pairs);
    assert_eq!(base.roots, shifted.roots);

    let base_elapsed = Array1::from_iter(base.mergers.iter().map(|merger| children[0].time - merger.time));
    let shifted_elapsed = Array1::from_iter(
      shifted
        .mergers
        .iter()
        .map(|merger| shifted_children[0].time - merger.time),
    );
    pretty_assert_abs_diff_eq!(base_elapsed, shifted_elapsed, epsilon = 1e-10);
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_sweep_plan_invariants(seed in any::<u64>()) {
      let children = vec![
        lineage(10.0, 0),
        lineage(9.5, 2),
        lineage(8.0, 1),
        lineage(7.0, 0),
        lineage(6.0, 3),
        lineage(2.0, 0),
      ];
      let t_stop = -5.0;
      let mut rng = get_random_number_generator(Some(seed));
      let plan = simulate_subtree(&children, t_stop, 0.7, &const_merger_rate(0.3), &mut rng).unwrap();

      assert_every_child_placed_once(&plan, children.len());

      // Time of every lineage the plan can name.
      let mut times: Vec<f64> = children.iter().map(|child| child.time).collect();
      times.extend(plan.mergers.iter().map(|merger| merger.time));

      let all_after_parent = plan.mergers.iter().all(|merger| merger.time > t_stop);
      prop_assert!(all_after_parent, "seed {seed}: merger at or past the parent bound {t_stop}: {:?}", plan.mergers);

      let all_references_valid = plan.mergers.iter().enumerate().all(|(index, merger)| {
        [merger.left, merger.right].into_iter().all(|lineage_id| {
          lineage_id < children.len() + index && merger.time <= times[lineage_id]
        })
      });
      prop_assert!(all_references_valid, "seed {seed}: merger references are invalid: {:?}", plan.mergers);

      // Sweep order requires merger times to decrease toward the parent.
      prop_assert!(
        plan.mergers.is_sorted_by(|first, second| first.time >= second.time),
        "seed {seed}: merger times are not monotone: {:?}",
        plan.mergers.iter().map(|m| m.time).collect::<Vec<_>>()
      );
    }

    #[test]
    fn test_prop_sweep_preserves_mutated_lineages_when_mutation_rate_is_zero(seed in any::<u64>()) {
      // With no mutation events, a mutated child is permanently ineligible to coalesce.
      let children = vec![
        lineage(10.0, 0),
        lineage(9.0, 0),
        lineage(8.0, 4),
        lineage(7.0, 0),
        lineage(6.0, 1),
      ];
      let mut rng = get_random_number_generator(Some(seed));
      let plan = simulate_subtree(&children, -50.0, 0.0, &const_merger_rate(2.0), &mut rng).unwrap();

      let mutated_children_survive = [2_usize, 4].into_iter().all(|child| plan.roots.contains(&child));
      prop_assert!(mutated_children_survive, "seed {seed}: a mutated child merged: {:?}", plan.mergers);
      assert_every_child_placed_once(&plan, children.len());
    }

    #[test]
    fn test_prop_sweep_resolves_fully_with_generous_window(seed in any::<u64>()) {
      let n_children = 8;
      let children: Vec<Lineage> = (0..n_children).map(|i| lineage(10.0 - i as f64, 0)).collect();
      let mut rng = get_random_number_generator(Some(seed));
      let plan = simulate_subtree(&children, -1.0e6, 0.0, &const_merger_rate(1.0), &mut rng).unwrap();

      prop_assert_eq!(plan.mergers.len(), n_children - 2);
      prop_assert_eq!(plan.roots.len(), 2);
    }
  }

  #[test]
  fn test_sweep_integrates_from_lineage_arrival_boundary() -> Result<(), Report> {
    let children = [lineage(10.0, 0), lineage(5.0, 0), lineage(0.0, 0)];
    let mut rng = get_random_number_generator(Some(3));

    let plan = simulate_subtree(&children, -100.0, 0.0, &const_merger_rate(1e6), &mut rng)?;

    let merger = plan.mergers.first().expect("two ready lineages must merge");
    assert!(
      4.99 < merger.time && merger.time < 5.0,
      "merger at {} must occur immediately after the second lineage arrives",
      merger.time
    );
    Ok(())
  }

  #[test]
  fn test_sweep_integrates_across_merger_rate_boundary() -> Result<(), Report> {
    let children = [lineage(10.0, 0), lineage(10.0, 0), lineage(10.0, 0)];
    // No merger hazard above year 5; a high hazard starts immediately below it.
    let merger_rate = PiecewiseConstantFn::new(array![5.0], array![1e6, 0.0]);
    let mut rng = get_random_number_generator(Some(3));

    let plan = simulate_subtree(&children, 0.0, 0.0, &merger_rate, &mut rng)?;

    let merger = plan
      .mergers
      .first()
      .expect("high post-boundary hazard must merge a pair");
    assert!(
      4.99 < merger.time && merger.time < 5.0,
      "merger at {} must occur immediately after crossing below year 5",
      merger.time
    );
    Ok(())
  }

  #[test]
  fn test_sweep_terminates_when_no_event_is_possible() -> Result<(), Report> {
    // No mutation rate and no merger rate: nothing can happen and nothing is pending once all
    // children are live. The sweep must return rather than spin on a zero-rate draw.
    let children = [lineage(3.0, 0), lineage(2.0, 1), lineage(1.0, 0)];
    let mut rng = get_random_number_generator(Some(1));

    let plan = simulate_subtree(&children, 0.0, 0.0, &const_merger_rate(0.0), &mut rng)?;

    assert!(plan.mergers.is_empty());
    assert_eq!(plan.roots.len(), 3);
    Ok(())
  }

  #[test]
  fn test_sweep_first_merger_waiting_time_matches_the_coalescent_rate() -> Result<(), Report> {
    // All children share a node time, so all are live from the start and the first merger
    // waits Exp((k - 1) * kappa).
    let k = 5;
    let kappa = 0.25;
    let expected_mean = 1.0 / ((k - 1) as f64 * kappa);

    #[allow(clippy::map_with_unused_argument_over_ranges)]
    let children: Vec<Lineage> = (0..k).map(|_| lineage(0.0, 0)).collect();
    let replicates = 20_000;

    let mut total = 0.0;
    for seed in 0..replicates {
      let mut rng = get_random_number_generator(Some(seed));
      let plan = simulate_subtree(&children, -1.0e9, 0.0, &const_merger_rate(kappa), &mut rng)?;
      let first = plan.mergers.first().expect("the window is effectively unbounded");
      total += -first.time;
    }
    let mean = total / replicates as f64;

    // Standard error of the mean is expected_mean / sqrt(replicates), about 0.7%.
    let relative_error = (mean - expected_mean).abs() / expected_mean;
    assert!(
      relative_error < 0.05,
      "mean first-merger wait {mean} deviates from the expected {expected_mean} by {:.1}%",
      relative_error * 100.0
    );
    Ok(())
  }
}
