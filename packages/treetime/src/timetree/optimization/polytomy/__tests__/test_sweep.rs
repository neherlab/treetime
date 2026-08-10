#[cfg(test)]
mod tests {
  use crate::timetree::optimization::polytomy::sweep::{Lineage, SubtreePlan, simulate_subtree};
  use eyre::Report;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use treetime_utils::assert_error;
  use treetime_utils::sync::random::get_random_number_generator;

  /// A merger rate that ignores time, so tests can reason about the sweep alone.
  fn const_merger_rate(rate: f64) -> impl Fn(f64) -> Result<f64, Report> {
    move |_time| Ok(rate)
  }

  /// A constant merger rate that also records the calendar time of every query, which is the
  /// sweep's position.
  fn recording_merger_rate(rate: f64, queries: &Mutex<Vec<f64>>) -> impl Fn(f64) -> Result<f64, Report> {
    move |time| {
      queries.lock().push(time);
      Ok(rate)
    }
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
  #[case::negative(-1.0,          "Polytomy merger rate must be finite and non-negative at calendar time 3.000000e0, got -1.000000e0")]
  #[case::nan(     f64::NAN,      "Polytomy merger rate must be finite and non-negative at calendar time 3.000000e0, got NaN")]
  #[case::infinite(f64::INFINITY, "Polytomy merger rate must be finite and non-negative at calendar time 3.000000e0, got inf")]
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
  fn test_sweep_plan_invariants_hold_across_seeds() -> Result<(), Report> {
    let children: Vec<Lineage> = vec![
      lineage(10.0, 0),
      lineage(9.5, 2),
      lineage(8.0, 1),
      lineage(7.0, 0),
      lineage(6.0, 3),
      lineage(2.0, 0),
    ];
    let t_stop = -5.0;

    for seed in 0..200 {
      let mut rng = get_random_number_generator(Some(seed));
      let plan = simulate_subtree(&children, t_stop, 0.7, &const_merger_rate(0.3), &mut rng)?;

      assert_every_child_placed_once(&plan, children.len());

      // Time of every lineage the plan can name.
      let mut times: Vec<f64> = children.iter().map(|child| child.time).collect();
      times.extend(plan.mergers.iter().map(|merger| merger.time));

      for (index, merger) in plan.mergers.iter().enumerate() {
        assert!(
          merger.time > t_stop,
          "seed {seed}: merger at {} is at or past the parent bound {t_stop}",
          merger.time
        );

        for lineage_id in [merger.left, merger.right] {
          assert!(
            lineage_id < children.len() + index,
            "seed {seed}: merger {index} references lineage {lineage_id}, which does not exist yet"
          );
          assert!(
            merger.time <= times[lineage_id],
            "seed {seed}: merger at {} is more recent than its child at {}",
            merger.time,
            times[lineage_id]
          );
        }
      }

      // Mergers are emitted oldest last.
      assert!(
        plan.mergers.is_sorted_by(|first, second| first.time >= second.time),
        "seed {seed}: merger times are not monotone: {:?}",
        plan.mergers.iter().map(|m| m.time).collect::<Vec<_>>()
      );
    }
    Ok(())
  }

  #[test]
  fn test_sweep_never_merges_a_lineage_that_still_has_mutations() -> Result<(), Report> {
    // With a zero mutation rate no substitution can ever be placed, so any child carrying one
    // is permanently ineligible to coalesce and must survive as a root.
    let children = vec![
      lineage(10.0, 0),
      lineage(9.0, 0),
      lineage(8.0, 4),
      lineage(7.0, 0),
      lineage(6.0, 1),
    ];

    for seed in 0..100 {
      let mut rng = get_random_number_generator(Some(seed));
      let plan = simulate_subtree(&children, -50.0, 0.0, &const_merger_rate(2.0), &mut rng)?;

      for mutated in [2_usize, 4] {
        assert!(
          plan.roots.contains(&mutated),
          "seed {seed}: child {mutated} carries mutations that can never be placed, so it cannot merge"
        );
      }
      assert_every_child_placed_once(&plan, children.len());
    }
    Ok(())
  }

  #[test]
  fn test_sweep_resumes_at_the_arrival_time_rather_than_the_drawn_time() -> Result<(), Report> {
    // A merger rate small enough that every draw overshoots the remaining window, so the sweep
    // is driven entirely by arrivals. The recorded query times are the sweep's positions.
    //
    // v0 resumes at the drawn time, skipping the interval between the arrival and that time.
    // Correct behaviour resumes at the arrival, so every query lands on a child's time.
    let children = [lineage(10.0, 0), lineage(9.0, 0), lineage(0.0, 0)];
    let queries = Mutex::new(Vec::new());
    let mut rng = get_random_number_generator(Some(3));

    let plan = simulate_subtree(&children, -100.0, 0.0, &recording_merger_rate(1e-9, &queries), &mut rng)?;

    let queries = queries.lock().clone();
    assert_eq!(
      queries,
      vec![10.0, 9.0, 0.0],
      "the sweep must step from arrival to arrival, not to the overshooting draw"
    );
    assert!(plan.mergers.is_empty(), "the rate is too low for a merger to fit");
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
  fn test_sweep_resolves_fully_when_the_window_is_generous() -> Result<(), Report> {
    // k mutation-free children and ample time: the sweep should run until only two lineages
    // remain, which for k children takes k - 2 mergers.
    let k = 8;
    let children: Vec<Lineage> = (0..k).map(|i| lineage(10.0 - i as f64, 0)).collect();

    for seed in 0..50 {
      let mut rng = get_random_number_generator(Some(seed));
      let plan = simulate_subtree(&children, -1.0e6, 0.0, &const_merger_rate(1.0), &mut rng)?;

      assert_eq!(plan.mergers.len(), k - 2, "seed {seed}: expected full resolution");
      assert_eq!(
        plan.roots.len(),
        2,
        "seed {seed}: a resolved polytomy leaves a bifurcation"
      );
    }
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
