#[cfg(test)]
mod tests {
  use crate::commands::mugration::input::MugrationInput;
  use crate::commands::mugration::run::{
    apply_pseudo_counts, compute_pi_from_weights, compute_pi_uniform, execute_mugration, validate_weight_coverage,
  };
  use crate::representation::discrete_states::DiscreteStates;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indexmap::IndexSet;
  use std::iter::once;
  use itertools::Itertools;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
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
    let graph = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };

    let input = MugrationInput {
      graph,
      traits,
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    };

    let result = execute_mugration(input)?;

    assert_eq!(o!("country"), result.gtr.attribute);
    assert_eq!(2, result.gtr.n_states);
    assert_eq!(vec_of_owned!["germany", "usa"], result.gtr.states);
    assert_eq!(2, result.gtr.pi.len());
    assert_abs_diff_eq!(1.0, result.gtr.pi.sum(), epsilon = 1e-10);
    assert!(
      result.gtr.mu > 0.1 && result.gtr.mu < 100.0,
      "mu should be in reasonable range for 2-state model: {}",
      result.gtr.mu
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
    let graph = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
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

    let input = MugrationInput {
      graph,
      traits,
      attribute: o!("country"),
      weights: Some(weights),
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    };

    let result = execute_mugration(input)?;

    assert_eq!(3, result.gtr.n_states);
    assert_eq!(vec_of_owned!["france", "germany", "usa"], result.gtr.states);
    assert_abs_diff_eq!(1.0, result.gtr.pi.sum(), epsilon = 1e-12);

    let total = 7.0;
    assert_abs_diff_eq!(1.0 / total, result.gtr.pi[0], epsilon = 1e-12); // france
    assert_abs_diff_eq!(4.0 / total, result.gtr.pi[1], epsilon = 1e-12); // germany
    assert_abs_diff_eq!(2.0 / total, result.gtr.pi[2], epsilon = 1e-12); // usa

    Ok(())
  }

  #[test]
  fn test_execute_mugration_with_weights_includes_unobserved_weight_states() -> Result<(), Report> {
    let graph = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };
    let weights = btreemap! {
      o!("usa") => 2.0,
      o!("germany") => 4.0,
      o!("france") => 1.0,
    };

    let input = MugrationInput {
      graph,
      traits,
      attribute: o!("country"),
      weights: Some(weights),
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    };

    let result = execute_mugration(input)?;

    assert_eq!(3, result.gtr.n_states);
    assert_eq!(vec_of_owned!["france", "germany", "usa"], result.gtr.states);
    assert_abs_diff_eq!(1.0, result.gtr.pi.sum(), epsilon = 1e-12);

    let total = 7.0;
    assert_abs_diff_eq!(1.0 / total, result.gtr.pi[0], epsilon = 1e-12); // france
    assert_abs_diff_eq!(4.0 / total, result.gtr.pi[1], epsilon = 1e-12); // germany
    assert_abs_diff_eq!(2.0 / total, result.gtr.pi[2], epsilon = 1e-12); // usa

    Ok(())
  }

  #[test]
  fn test_execute_mugration_with_pseudo_counts() -> Result<(), Report> {
    let graph = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
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

    let input = MugrationInput {
      graph,
      traits,
      attribute: o!("country"),
      weights: Some(weights),
      missing_data: o!("?"),
      pc: Some(1.0),
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    };

    let result = execute_mugration(input)?;

    // With iterative GTR and fixed_pi from weights, final pi reflects weight proportions:
    // weights = {usa: 3, germany: 1, france: 1}, normalized = [0.2, 0.2, 0.6]
    // Pseudo-counts affect the initial model but fixed_pi keeps weight proportions in refinement.
    assert_abs_diff_eq!(1.0, result.gtr.pi.sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(0.2, result.gtr.pi[0], epsilon = 1e-10); // france
    assert_abs_diff_eq!(0.2, result.gtr.pi[1], epsilon = 1e-10); // germany
    assert_abs_diff_eq!(0.6, result.gtr.pi[2], epsilon = 1e-10); // usa

    Ok(())
  }

  #[test]
  fn test_execute_mugration_sampling_bias_correction() -> Result<(), Report> {
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };

    let base_input = MugrationInput {
      graph: nwk_read_str("(A:0.1,B:0.2)root;")?,
      traits: traits.clone(),
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    };
    let base_result = execute_mugration(base_input)?;
    let base_mu = base_result.gtr.mu;

    let corrected_input = MugrationInput {
      graph: nwk_read_str("(A:0.1,B:0.2)root;")?,
      traits,
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: Some(2.0),
    };
    let corrected_result = execute_mugration(corrected_input)?;

    // sampling_bias_correction multiplies mu after iterative refinement
    assert_abs_diff_eq!(corrected_result.gtr.mu, base_mu * 2.0, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_execute_mugration_rejects_single_state() {
    let graph = nwk_read_str("(A:0.1,B:0.2)root;").unwrap();
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("usa"),
    };

    let input = MugrationInput {
      graph,
      traits,
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    };

    let result = execute_mugration(input);
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

    let result_no_iter = execute_mugration(MugrationInput {
      graph: nwk_read_str(tree)?,
      traits: traits.clone(),
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 0,
      sampling_bias_correction: None,
    })?;

    let result_with_iter = execute_mugration(MugrationInput {
      graph: nwk_read_str(tree)?,
      traits,
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    })?;

    // Iterative refinement must change the model: mu and/or pi should differ
    let mu_changed = (result_no_iter.gtr.mu - result_with_iter.gtr.mu).abs() > 1e-6;
    let pi_changed = result_no_iter
      .gtr
      .pi
      .iter()
      .zip(result_with_iter.gtr.pi.iter())
      .any(|(a, b)| (a - b).abs() > 1e-6);
    assert!(
      mu_changed || pi_changed,
      "iterative refinement must change the model: mu_0={}, mu_5={}, pi_0={:?}, pi_5={:?}",
      result_no_iter.gtr.mu,
      result_with_iter.gtr.mu,
      result_no_iter.gtr.pi,
      result_with_iter.gtr.pi
    );

    // Both models must have valid parameters
    assert!(result_with_iter.gtr.mu > 0.0);
    assert_abs_diff_eq!(result_with_iter.gtr.pi.sum(), 1.0, epsilon = 1e-10);
    assert!(result_with_iter.gtr.pi.iter().all(|&p| p > 0.0));

    Ok(())
  }

  #[test]
  fn test_iterative_refinement_pi_reflects_data() -> Result<(), Report> {
    // 4 tips: 3 usa, 1 germany. Without fixed_pi, iterative refinement
    // should shift pi toward observed frequencies (higher weight for usa).
    let tree = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;";
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("usa"),
      o!("C") => o!("usa"),
      o!("D") => o!("germany"),
    };

    let result = execute_mugration(MugrationInput {
      graph: nwk_read_str(tree)?,
      traits,
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
    })?;

    // pi for usa (index 1, alphabetically: germany=0, usa=1) should be > 0.5
    let pi_usa = result.gtr.pi[1];
    assert!(
      pi_usa > 0.5,
      "pi[usa] should reflect 3/4 observed frequency, got {pi_usa:.4}"
    );

    Ok(())
  }

  #[test]
  fn test_zero_iterations_preserves_initial_model() -> Result<(), Report> {
    let result = execute_mugration(MugrationInput {
      graph: nwk_read_str("(A:0.1,B:0.2)root;")?,
      traits: btreemap! { o!("A") => o!("usa"), o!("B") => o!("germany") },
      attribute: o!("country"),
      weights: None,
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 0,
      sampling_bias_correction: None,
    })?;

    // With 0 iterations and no weights, initial model has uniform pi
    // The infer+optimize cycle still runs once (initial inference), but with
    // iterations=0 the loop body doesn't execute.
    assert_abs_diff_eq!(result.gtr.pi.sum(), 1.0, epsilon = 1e-10);
    assert!(result.gtr.mu > 0.0);

    Ok(())
  }
}
