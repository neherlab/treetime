#[cfg(test)]
mod tests {
  use crate::alphabet::state_alphabet::StateAlphabet;
  use crate::mugration::pipeline::{
    MugrationParams, apply_pseudo_counts, compute_pi_from_weights, compute_pi_uniform, validate_weight_coverage,
  };
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use helpers::{confidence_of, mugration_params, run_mugration_case, trait_assignments_by_name};
  use indexmap::IndexSet;
  use itertools::Itertools;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::iter::once;
  use treetime_utils::{assert_error, o, vec_of_owned};

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
  fn test_run_compute_pi_from_weights_normalizes() -> Result<(), Report> {
    let states = StateAlphabet::new(["usa", "germany"])?;
    let weights = btreemap! {
      o!("usa") => 3.0,
      o!("germany") => 1.0,
    };

    let pi = compute_pi_from_weights(&states, &weights);

    assert_eq!(2, pi.len());
    assert_abs_diff_eq!(1.0, pi.sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(0.25, pi[0], epsilon = 1e-10);
    assert_abs_diff_eq!(0.75, pi[1], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_run_compute_pi_from_weights_uses_mean_for_missing() -> Result<(), Report> {
    let states = StateAlphabet::new(["usa", "germany", "france"])?;
    let weights = btreemap! {
      o!("usa") => 2.0,
      o!("germany") => 4.0,
    };

    let pi = compute_pi_from_weights(&states, &weights);

    assert_eq!(3, pi.len());
    assert_abs_diff_eq!(1.0, pi.sum(), epsilon = 1e-10);
    let total = 9.0;
    assert_abs_diff_eq!(3.0 / total, pi[0], epsilon = 1e-10);
    assert_abs_diff_eq!(4.0 / total, pi[1], epsilon = 1e-10);
    assert_abs_diff_eq!(2.0 / total, pi[2], epsilon = 1e-10);
    Ok(())
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
  fn test_run_apply_pseudo_counts_adds_count_and_renormalizes() {
    let pi = array![0.25, 0.75];
    let pc = 0.5;

    let result = apply_pseudo_counts(pi, pc);

    assert_eq!(2, result.len());
    assert_abs_diff_eq!(1.0, result.sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(0.375, result[0], epsilon = 1e-10);
    assert_abs_diff_eq!(0.625, result[1], epsilon = 1e-10);
  }

  #[test]
  fn test_run_apply_pseudo_counts_preserves_normalization() {
    let pi = array![0.1, 0.2, 0.3, 0.4];
    let pc = 1.0;

    let result = apply_pseudo_counts(pi, pc);

    assert_abs_diff_eq!(1.0, result.sum(), epsilon = 1e-10);
  }

  #[test]
  fn test_run_mugration_simple_tree() -> Result<(), Report> {
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };
    let (output, names) = run_mugration_case("(A:0.1,B:0.2)root;", &traits, None, &mugration_params())?;

    assert_eq!(2, output.states.len());
    assert_eq!(
      vec_of_owned!["germany", "usa"],
      output.states.iter().map(|s| s.to_owned()).collect_vec()
    );
    assert_eq!(2, output.gtr.pi().len());
    assert_abs_diff_eq!(1.0, output.gtr.pi().sum(), epsilon = 1e-10);
    assert!(
      output.gtr.mu() > 0.1 && output.gtr.mu() < 100.0,
      "mu should be in reasonable range for 2-state model: {}",
      output.gtr.mu()
    );

    let assignments = trait_assignments_by_name(&output, &names);
    assert_eq!(3, assignments.len());
    assert_eq!(Some(&o!("usa")), assignments.get("root"));
    assert_eq!(Some(&o!("usa")), assignments.get("A"));
    assert_eq!(Some(&o!("germany")), assignments.get("B"));

    let root_confidence = confidence_of(&output, &names, "root");
    assert_abs_diff_eq!(1.0, root_confidence.sum(), epsilon = 1e-12);

    Ok(())
  }

  #[test]
  fn test_run_mugration_with_weights() -> Result<(), Report> {
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
    let (output, _names) = run_mugration_case("(A:0.1,B:0.2,C:0.3)root;", &traits, Some(weights), &mugration_params())?;

    assert_eq!(3, output.states.len());
    assert_eq!(
      vec_of_owned!["france", "germany", "usa"],
      output.states.iter().map(|s| s.to_owned()).collect_vec()
    );
    assert_abs_diff_eq!(1.0, output.gtr.pi().sum(), epsilon = 1e-12);

    let total = 7.0;
    assert_abs_diff_eq!(1.0 / total, output.gtr.pi()[0], epsilon = 1e-12);
    assert_abs_diff_eq!(4.0 / total, output.gtr.pi()[1], epsilon = 1e-12);
    assert_abs_diff_eq!(2.0 / total, output.gtr.pi()[2], epsilon = 1e-12);

    Ok(())
  }

  #[test]
  fn test_run_mugration_with_weights_includes_unobserved_weight_states() -> Result<(), Report> {
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };
    let weights = btreemap! {
      o!("usa") => 2.0,
      o!("germany") => 4.0,
      o!("france") => 1.0,
    };
    let (output, _names) = run_mugration_case("(A:0.1,B:0.2)root;", &traits, Some(weights), &mugration_params())?;

    assert_eq!(3, output.states.len());
    assert_eq!(
      vec_of_owned!["france", "germany", "usa"],
      output.states.iter().map(|s| s.to_owned()).collect_vec()
    );
    assert_abs_diff_eq!(1.0, output.gtr.pi().sum(), epsilon = 1e-12);

    let total = 7.0;
    assert_abs_diff_eq!(1.0 / total, output.gtr.pi()[0], epsilon = 1e-12);
    assert_abs_diff_eq!(4.0 / total, output.gtr.pi()[1], epsilon = 1e-12);
    assert_abs_diff_eq!(2.0 / total, output.gtr.pi()[2], epsilon = 1e-12);

    Ok(())
  }

  #[test]
  fn test_run_mugration_with_pseudo_counts() -> Result<(), Report> {
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
    let (output, _names) = run_mugration_case(
      "(A:0.1,B:0.2,C:0.3)root;",
      &traits,
      Some(weights),
      &MugrationParams {
        pc: Some(1.0),
        ..mugration_params()
      },
    )?;

    assert_abs_diff_eq!(1.0, output.gtr.pi().sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(0.2, output.gtr.pi()[0], epsilon = 1e-10);
    assert_abs_diff_eq!(0.2, output.gtr.pi()[1], epsilon = 1e-10);
    assert_abs_diff_eq!(0.6, output.gtr.pi()[2], epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_run_mugration_smooth_initial_pi_preserves_fixed_equilibrium() -> Result<(), Report> {
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
    let (output, _names) = run_mugration_case(
      "(A:0.1,B:0.2,C:0.3)root;",
      &traits,
      Some(weights),
      &MugrationParams {
        pc: Some(1.0),
        smooth_initial_pi: true,
        ..mugration_params()
      },
    )?;

    let total = 5.0;
    assert_abs_diff_eq!(1.0, output.gtr.pi().sum(), epsilon = 1e-10);
    assert_abs_diff_eq!(1.0 / total, output.gtr.pi()[0], epsilon = 1e-10);
    assert_abs_diff_eq!(1.0 / total, output.gtr.pi()[1], epsilon = 1e-10);
    assert_abs_diff_eq!(3.0 / total, output.gtr.pi()[2], epsilon = 1e-10);

    Ok(())
  }

  #[test]
  fn test_run_mugration_filter_uninformative_root_changes_inferred_equilibrium() -> Result<(), Report> {
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };

    let (v0, _v0_names) = run_mugration_case(
      "(A:0.1,B:0.1)root;",
      &traits,
      None,
      &MugrationParams {
        iterations: 0,
        ..mugration_params()
      },
    )?;
    let (filtered, _filtered_names) = run_mugration_case(
      "(A:0.1,B:0.1)root;",
      &traits,
      None,
      &MugrationParams {
        iterations: 0,
        filter_uninformative_root: true,
        ..mugration_params()
      },
    )?;

    let pi_v0 = &v0.gtr.pi();
    let pi_filtered = &filtered.gtr.pi();
    assert!(
      (pi_v0[0] - pi_filtered[0]).abs() > 1e-9,
      "uninformative-root filtering must change the inferred equilibrium: v0 pi={pi_v0:?}, filtered pi={pi_filtered:?}"
    );

    Ok(())
  }

  #[test]
  fn test_run_mugration_sampling_bias_correction() -> Result<(), Report> {
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
    };

    let (base, _base_names) = run_mugration_case("(A:0.1,B:0.2)root;", &traits, None, &mugration_params())?;
    let base_mu = base.gtr.mu();

    let (corrected, _corrected_names) = run_mugration_case(
      "(A:0.1,B:0.2)root;",
      &traits,
      None,
      &MugrationParams {
        sampling_bias_correction: Some(2.0),
        ..mugration_params()
      },
    )?;

    assert_abs_diff_eq!(corrected.gtr.mu(), base_mu * 2.0, epsilon = 1e-6);

    Ok(())
  }

  #[test]
  fn test_run_mugration_rejects_single_state() {
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("usa"),
    };

    let result = run_mugration_case("(A:0.1,B:0.2)root;", &traits, None, &mugration_params());
    assert_error!(
      result,
      "A discrete trait needs at least 2 different states, found 1 state: usa. \
       Check that the trait column holds at least 2 different known values."
    );
  }

  #[test]
  fn test_run_mugration_imputes_leaf_missing_from_metadata() -> Result<(), Report> {
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("C") => o!("germany"),
    };
    let (output, names) = run_mugration_case("(A:0.1,B:0.2,C:0.3)root;", &traits, None, &mugration_params())?;

    let assignments = trait_assignments_by_name(&output, &names);
    assert_eq!(
      vec_of_owned!["A", "B", "C", "root"],
      assignments.keys().cloned().sorted().collect_vec()
    );
    let leaf_b_confidence = confidence_of(&output, &names, "B");
    assert_abs_diff_eq!(1.0, leaf_b_confidence.sum(), epsilon = 1e-12);

    Ok(())
  }

  #[test]
  fn test_iterative_refinement_changes_model() -> Result<(), Report> {
    let tree = "(A:0.1,(B:0.05,C:0.15)BC:0.2)root;";
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("usa"),
    };

    let (no_iter, _no_iter_names) = run_mugration_case(
      tree,
      &traits,
      None,
      &MugrationParams {
        iterations: 0,
        ..mugration_params()
      },
    )?;
    let (with_iter, _with_iter_names) = run_mugration_case(tree, &traits, None, &mugration_params())?;

    let mu_changed = (no_iter.gtr.mu() - with_iter.gtr.mu()).abs() > 1e-6;
    let pi_changed = no_iter
      .gtr
      .pi()
      .iter()
      .zip(with_iter.gtr.pi().iter())
      .any(|(a, b)| (a - b).abs() > 1e-6);
    assert!(
      mu_changed || pi_changed,
      "iterative refinement must change the model: mu_0={}, mu_5={}, pi_0={:?}, pi_5={:?}",
      no_iter.gtr.mu(),
      with_iter.gtr.mu(),
      no_iter.gtr.pi(),
      with_iter.gtr.pi()
    );

    assert!(with_iter.gtr.mu() > 0.0);
    assert_abs_diff_eq!(with_iter.gtr.pi().sum(), 1.0, epsilon = 1e-10);
    assert!(with_iter.gtr.pi().iter().all(|&p| p > 0.0));

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

    let (output, _names) = run_mugration_case(tree, &traits, None, &mugration_params())?;

    let pi_usa = output.gtr.pi()[1];
    assert!(
      pi_usa > 0.5,
      "pi[usa] should reflect 3/4 observed frequency, got {pi_usa:.4}"
    );

    Ok(())
  }

  #[test]
  fn test_zero_iterations_preserves_initial_model() -> Result<(), Report> {
    let traits = btreemap! { o!("A") => o!("usa"), o!("B") => o!("germany") };
    let (output, _names) = run_mugration_case(
      "(A:0.1,B:0.2)root;",
      &traits,
      None,
      &MugrationParams {
        iterations: 0,
        ..mugration_params()
      },
    )?;

    assert_abs_diff_eq!(output.gtr.pi().sum(), 1.0, epsilon = 1e-10);
    assert!(output.gtr.mu() > 0.0);

    Ok(())
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::mugration::pipeline::{MugrationInput, MugrationOutput, MugrationParams, run};
    use crate::mugration::trait_pairing::{TraitPairing, pair_trait_rows};
    use crate::progress::NoopProgress;
    use eyre::Report;
    use indexmap::IndexMap;
    use ndarray::Array1;
    use std::collections::BTreeMap;
    use treetime_graph::node_names::NodeNames;

    use treetime_io::nwk::nwk_read;
    use treetime_utils::o;

    pub(super) fn mugration_params() -> MugrationParams {
      MugrationParams {
        missing_data: o!("?"),
        pc: None,
        missing_weights_threshold: 0.5,
        iterations: 5,
        sampling_bias_correction: None,
        smooth_initial_pi: false,
        filter_uninformative_root: false,
      }
    }

    pub(super) fn run_mugration_case(
      nwk: &str,
      traits: &BTreeMap<String, String>,
      weights: Option<BTreeMap<String, f64>>,
      params: &MugrationParams,
    ) -> Result<(MugrationOutput, NodeNames), Report> {
      let nwk_parsed = nwk_read(nwk.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let TraitPairing {
        traits,
        observed_values,
        ..
      } = pair_trait_rows(traits.clone(), &graph, &names);
      let input = MugrationInput {
        graph,
        traits,
        traits_source: o!("metadata"),
        observed_values,
        weights,
        branch_lengths,
      };
      let output = run(params, input, &names, &NoopCancel, &NoopProgress).map_err(|err| err.into_report())?;
      Ok((output, names))
    }

    pub(super) fn trait_assignments_by_name(output: &MugrationOutput, names: &NodeNames) -> IndexMap<String, String> {
      output
        .graph
        .get_nodes()
        .map(|node| {
          let key = node.key();
          let name = names[key].to_owned();
          let state = output.states.name(output.posteriors[&key].state()).unwrap();
          (name, state.to_owned())
        })
        .collect()
    }

    pub(super) fn confidence_of(output: &MugrationOutput, names: &NodeNames, node_name: &str) -> Array1<f64> {
      let key = output
        .graph
        .get_nodes()
        .map(|node| node.key())
        .find(|key| &names[*key] == node_name)
        .unwrap_or_else(|| panic!("missing node '{node_name}'"));
      output.posteriors[&key].profile().clone()
    }
  }
}
