#[cfg(test)]
mod tests {
  use treetime::mugration::mugration::{
    apply_pseudo_counts, compute_pi_from_weights, compute_pi_uniform, execute_mugration, validate_weight_coverage,
  };
  use treetime::partition::storage::discrete::DiscreteStates;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indexmap::IndexSet;
  use itertools::Itertools;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::iter::once;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_run_validate_weight_coverage_rejects_above_threshold() {
    let unique_values: IndexSet<String> = [o!("usa"), o!("germany"), o!("france"), o!("italy")]
      .into_iter()
      .collect();
    let weights_keys: IndexSet<String> = once(o!("usa")).collect();
    let missing_data = "?";
    let threshold = 0.5;

    let result = validate_weight_coverage(&unique_values, &weights_keys, missing_data, threshold);
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("too many discrete attributes missing"));
    assert!(err_msg.contains("0.75"));
  }

  #[test]
  fn test_run_validate_weight_coverage_accepts_at_threshold() {
    let unique_values: IndexSet<String> = [o!("usa"), o!("germany")].into_iter().collect();
    let weights_keys: IndexSet<String> = once(o!("usa")).collect();
    let missing_data = "?";
    let threshold = 0.5;

    let coverage = validate_weight_coverage(&unique_values, &weights_keys, missing_data, threshold).unwrap();
    let expected_missing: IndexSet<String> = once(o!("germany")).collect();
    assert_eq!(expected_missing, coverage.missing_values);
    assert_abs_diff_eq!(0.5, coverage.missing_ratio, epsilon = 1e-10);
  }

  #[test]
  fn test_run_validate_weight_coverage_excludes_missing_data_marker() {
    let unique_values: IndexSet<String> = [o!("usa"), o!("?")].into_iter().collect();
    let weights_keys: IndexSet<String> = once(o!("usa")).collect();
    let missing_data = "?";
    let threshold = 0.5;

    let coverage = validate_weight_coverage(&unique_values, &weights_keys, missing_data, threshold).unwrap();
    assert!(coverage.missing_values.is_empty());
    assert_abs_diff_eq!(0.0, coverage.missing_ratio, epsilon = 1e-10);
  }

  #[test]
  fn test_run_validate_weight_coverage_full_coverage() {
    let unique_values: IndexSet<String> = [o!("usa"), o!("germany")].into_iter().collect();
    let weights_keys: IndexSet<String> = [o!("usa"), o!("germany")].into_iter().collect();
    let missing_data = "?";
    let threshold = 0.5;

    let coverage = validate_weight_coverage(&unique_values, &weights_keys, missing_data, threshold).unwrap();
    assert!(coverage.missing_values.is_empty());
    assert_abs_diff_eq!(0.0, coverage.missing_ratio, epsilon = 1e-10);
  }

  #[test]
  fn test_run_compute_pi_from_weights_normalizes() {
    let states = DiscreteStates::from_values(["usa", "germany"].into_iter(), "?");
    let weights = btreemap! {
      o!("usa") => 3.0,
      o!("germany") => 1.0,
    };

    let pi = compute_pi_from_weights(&states, &weights);

    assert_eq!(2, pi.len());
    assert_abs_diff_eq!(1.0, pi.sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(0.25, pi[0], epsilon = 1e-10);
    assert_abs_diff_eq!(0.75, pi[1], epsilon = 1e-10);
  }

  #[test]
  fn test_run_compute_pi_from_weights_uses_mean_for_missing() {
    // DiscreteStates sorts alphabetically: france, germany, usa
    let states = DiscreteStates::from_values(["usa", "germany", "france"].into_iter(), "?");
    let weights = btreemap! {
      o!("usa") => 2.0,
      o!("germany") => 4.0,
    };

    let pi = compute_pi_from_weights(&states, &weights);

    assert_eq!(3, pi.len());
    assert_abs_diff_eq!(1.0, pi.sum(), epsilon = 1e-10);
    // mean_weight = (2.0 + 4.0) / 2 = 3.0
    // total = france(3.0) + germany(4.0) + usa(2.0) = 9.0
    let total = 9.0;
    assert_abs_diff_eq!(3.0 / total, pi[0], epsilon = 1e-10); // france (mean)
    assert_abs_diff_eq!(4.0 / total, pi[1], epsilon = 1e-10); // germany
    assert_abs_diff_eq!(2.0 / total, pi[2], epsilon = 1e-10); // usa
  }

  #[test]
  fn test_run_compute_pi_uniform() {
    let pi = compute_pi_uniform(4);

    assert_eq!(4, pi.len());
    assert_abs_diff_eq!(1.0, pi.sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(0.25, pi[0], epsilon = 1e-10);
    assert_abs_diff_eq!(0.25, pi[1], epsilon = 1e-10);
    assert_abs_diff_eq!(0.25, pi[2], epsilon = 1e-10);
    assert_abs_diff_eq!(0.25, pi[3], epsilon = 1e-10);
  }

  #[test]
  fn test_run_apply_pseudo_counts_with_value() {
    let pi = array![0.25, 0.75];
    let pc = Some(0.5);

    let result = apply_pseudo_counts(pi, pc);

    assert_eq!(2, result.len());
    assert_abs_diff_eq!(1.0, result.sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(0.375, result[0], epsilon = 1e-10);
    assert_abs_diff_eq!(0.625, result[1], epsilon = 1e-10);
  }

  #[test]
  fn test_run_apply_pseudo_counts_without_value() {
    let pi = array![0.25, 0.75];
    let pc = None;

    let result = apply_pseudo_counts(pi.clone(), pc);

    assert_abs_diff_eq!(pi, result, epsilon = 1e-10);
  }

  #[test]
  fn test_run_apply_pseudo_counts_preserves_normalization() {
    let pi = array![0.1, 0.2, 0.3, 0.4];
    let pc = Some(1.0);

    let result = apply_pseudo_counts(pi, pc);

    assert_abs_diff_eq!(1.0, result.sum(), epsilon = 1e-10);
  }

  #[test]
  fn test_execute_mugration_simple_tree() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };

    let names_tt_14 = names;
    let (result, maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_14,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      5,
      None,
      false,
      false,
    )?;

    assert_eq!(o!("country"), result.traits.attribute);
    assert_eq!(2, maps.n_states);
    assert_eq!(
      vec_of_owned!["germany", "usa"],
      maps.states.iter().map(|s| s.to_owned()).collect_vec()
    );
    assert_eq!(2, maps.gtr.pi.len());
    assert_abs_diff_eq!(1.0, maps.gtr.pi.sum(), epsilon = 1e-10);
    assert!(
      maps.gtr.mu > 0.1 && maps.gtr.mu < 100.0,
      "mu should be in reasonable range for 2-state model: {}",
      maps.gtr.mu
    );

    assert_eq!(3, result.traits.assignments.len());
    assert_eq!(Some(&o!("usa")), result.traits.assignments.get("root"));
    assert_eq!(Some(&o!("usa")), result.traits.assignments.get("A"));
    assert_eq!(Some(&o!("germany")), result.traits.assignments.get("B"));

    assert_eq!(vec_of_owned!["germany", "usa"], result.confidence.states);
    assert_eq!(3, result.confidence.rows.len());

    let root_confidence = result.confidence.rows.iter().find(|r| r.node == "root").unwrap();
    assert_abs_diff_eq!(1.0, root_confidence.profile.sum(), epsilon = 1e-12);

    Ok(())
  }

  #[test]
  fn test_execute_mugration_with_weights() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("france"),
    };
    let weights = btreemap! {
      o!("usa") => 2.0,
      o!("germany") => 4.0,
      o!("france") => 1.0,
    };

    let names_tt_13 = names;
    let (_result, maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_13,
      &branch_lengths,
      &traits,
      "country",
      Some(&weights),
      "?",
      None,
      0.5,
      5,
      None,
      false,
      false,
    )?;

    assert_eq!(3, maps.n_states);
    assert_eq!(
      vec_of_owned!["france", "germany", "usa"],
      maps.states.iter().map(|s| s.to_owned()).collect_vec()
    );
    assert_abs_diff_eq!(1.0, maps.gtr.pi.sum(), epsilon = 1e-12);

    let total = 7.0;
    assert_abs_diff_eq!(1.0 / total, maps.gtr.pi[0], epsilon = 1e-12); // france
    assert_abs_diff_eq!(4.0 / total, maps.gtr.pi[1], epsilon = 1e-12); // germany
    assert_abs_diff_eq!(2.0 / total, maps.gtr.pi[2], epsilon = 1e-12); // usa

    Ok(())
  }

  #[test]
  fn test_execute_mugration_with_weights_includes_unobserved_weight_states() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };
    let weights = btreemap! {
      o!("usa") => 2.0,
      o!("germany") => 4.0,
      o!("france") => 1.0,
    };

    let names_tt_12 = names;
    let (_result, maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_12,
      &branch_lengths,
      &traits,
      "country",
      Some(&weights),
      "?",
      None,
      0.5,
      5,
      None,
      false,
      false,
    )?;

    assert_eq!(3, maps.n_states);
    assert_eq!(
      vec_of_owned!["france", "germany", "usa"],
      maps.states.iter().map(|s| s.to_owned()).collect_vec()
    );
    assert_abs_diff_eq!(1.0, maps.gtr.pi.sum(), epsilon = 1e-12);

    let total = 7.0;
    assert_abs_diff_eq!(1.0 / total, maps.gtr.pi[0], epsilon = 1e-12); // france
    assert_abs_diff_eq!(4.0 / total, maps.gtr.pi[1], epsilon = 1e-12); // germany
    assert_abs_diff_eq!(2.0 / total, maps.gtr.pi[2], epsilon = 1e-12); // usa

    Ok(())
  }

  #[test]
  fn test_execute_mugration_with_pseudo_counts() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("france"),
    };
    let weights = btreemap! {
      o!("usa") => 3.0,
      o!("germany") => 1.0,
      o!("france") => 1.0,
    };

    let names_tt_11 = names;
    let (_result, maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_11,
      &branch_lengths,
      &traits,
      "country",
      Some(&weights),
      "?",
      Some(1.0),
      0.5,
      5,
      None,
      false,
      false,
    )?;

    assert_abs_diff_eq!(1.0, maps.gtr.pi.sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(0.2, maps.gtr.pi[0], epsilon = 1e-10); // france
    assert_abs_diff_eq!(0.2, maps.gtr.pi[1], epsilon = 1e-10); // germany
    assert_abs_diff_eq!(0.6, maps.gtr.pi[2], epsilon = 1e-10); // usa

    Ok(())
  }

  #[test]
  fn test_execute_mugration_smooth_initial_pi_preserves_fixed_equilibrium() -> Result<(), Report> {
    // Smoothing flattens only the initial pi used for the first reconstruction
    // pass; the final equilibrium stays pinned to the raw weight-derived
    // fixed_pi, so the returned model matches the unsmoothed equilibrium.
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("france"),
    };
    let weights = btreemap! {
      o!("usa") => 3.0,
      o!("germany") => 1.0,
      o!("france") => 1.0,
    };

    let names_tt_10 = names;
    let (_result, maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_10,
      &branch_lengths,
      &traits,
      "country",
      Some(&weights),
      "?",
      Some(1.0),
      0.5,
      5,
      None,
      true,
      false,
    )?;

    let total = 5.0;
    assert_abs_diff_eq!(1.0, maps.gtr.pi.sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(1.0 / total, maps.gtr.pi[0], epsilon = 1e-10); // france
    assert_abs_diff_eq!(1.0 / total, maps.gtr.pi[1], epsilon = 1e-10); // germany
    assert_abs_diff_eq!(3.0 / total, maps.gtr.pi[2], epsilon = 1e-10); // usa

    Ok(())
  }

  #[test]
  fn test_execute_mugration_filter_uninformative_root_changes_inferred_equilibrium() -> Result<(), Report> {
    // Symmetric two-leaf tree: the root posterior is exactly uniform [0.5, 0.5].
    // With filtering off (v0) the root's argmax state is folded into the
    // equilibrium-frequency prior; with filtering on the uniform root is
    // skipped. With no weights (pi inferred from counts) the two policies must
    // therefore yield different equilibrium frequencies.
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };

    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let names_tt_9 = names;
    let (_v0, v0_maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_9,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      0,
      None,
      false,
      false,
    )?;
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let names_tt_8 = names;
    let (_filtered, filtered_maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_8,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      0,
      None,
      false,
      true,
    )?;

    let pi_v0 = &v0_maps.gtr.pi;
    let pi_filtered = &filtered_maps.gtr.pi;
    assert!(
      (pi_v0[0] - pi_filtered[0]).abs() > 1e-9,
      "uninformative-root filtering must change the inferred equilibrium: v0 pi={pi_v0:?}, filtered pi={pi_filtered:?}"
    );

    Ok(())
  }

  #[test]
  fn test_execute_mugration_sampling_bias_correction() -> Result<(), Report> {
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };

    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let names_tt_7 = names;
    let (_base_result, base_maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_7,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      5,
      None,
      false,
      false,
    )?;
    let base_mu = base_maps.gtr.mu;

    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let names_tt_6 = names;
    let (_corrected_result, corrected_maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_6,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      5,
      Some(2.0),
      false,
      false,
    )?;

    assert_abs_diff_eq!(corrected_maps.gtr.mu, base_mu * 2.0, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_execute_mugration_rejects_single_state() {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;").unwrap();
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("usa"),
    };

    let names_tt_5 = names;
    let result = execute_mugration(
      graph,
      &confidences,
      &names_tt_5,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      5,
      None,
      false,
      false,
    );
    assert!(result.is_err());
    let err_msg = result.unwrap_err().to_string();
    assert!(err_msg.contains("only 1 discrete attributes"));
    assert!(err_msg.contains("At least 2 are required"));
  }

  #[test]
  fn test_iterative_refinement_changes_model() -> Result<(), Report> {
    let tree = "(A:0.1,(B:0.05,C:0.15)BC:0.2)root;";
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("usa"),
    };

    let nwk_parsed = nwk_read_str(tree)?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let names_tt_4 = names;
    let (_result_no_iter, no_iter_maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_4,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      0,
      None,
      false,
      false,
    )?;

    let nwk_parsed = nwk_read_str(tree)?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let names_tt_3 = names;
    let (_result_with_iter, with_iter_maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_3,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      5,
      None,
      false,
      false,
    )?;

    let mu_changed = (no_iter_maps.gtr.mu - with_iter_maps.gtr.mu).abs() > 1e-6;
    let pi_changed = no_iter_maps
      .gtr
      .pi
      .iter()
      .zip(with_iter_maps.gtr.pi.iter())
      .any(|(a, b)| (a - b).abs() > 1e-6);
    assert!(
      mu_changed || pi_changed,
      "iterative refinement must change the model: mu_0={}, mu_5={}, pi_0={:?}, pi_5={:?}",
      no_iter_maps.gtr.mu,
      with_iter_maps.gtr.mu,
      no_iter_maps.gtr.pi,
      with_iter_maps.gtr.pi
    );

    assert!(with_iter_maps.gtr.mu > 0.0);
    assert_abs_diff_eq!(with_iter_maps.gtr.pi.sum(), 1.0, epsilon = 1e-10);
    assert!(with_iter_maps.gtr.pi.iter().all(|&p| p > 0.0));

    Ok(())
  }

  #[test]
  fn test_iterative_refinement_pi_reflects_data() -> Result<(), Report> {
    let tree = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;";
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("usa"),
      o!("C") => o!("usa"),
      o!("D") => o!("germany"),
    };

    let nwk_parsed = nwk_read_str(tree)?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let names_tt_2 = names;
    let (_result, maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_2,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      5,
      None,
      false,
      false,
    )?;

    let pi_usa = maps.gtr.pi[1];
    assert!(
      pi_usa > 0.5,
      "pi[usa] should reflect 3/4 observed frequency, got {pi_usa:.4}"
    );

    Ok(())
  }

  #[test]
  fn test_zero_iterations_preserves_initial_model() -> Result<(), Report> {
    let traits = btreemap! { o!("A") => o!("usa"), o!("B") => o!("germany") };
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let names_tt_1 = names;
    let (_result, maps) = execute_mugration(
      graph,
      &confidences,
      &names_tt_1,
      &branch_lengths,
      &traits,
      "country",
      None,
      "?",
      None,
      0.5,
      0,
      None,
      false,
      false,
    )?;

    assert_abs_diff_eq!(maps.gtr.pi.sum(), 1.0, epsilon = 1e-10);
    assert!(maps.gtr.mu > 0.0);

    Ok(())
  }
}
