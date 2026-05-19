#[cfg(test)]
mod tests {
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  use helpers::{load_gm_mugration_inputs, load_gm_mugration_outputs, run_gm_mugration_case};

  // Golden master tests for mugration discrete trait reconstruction.
  //
  // Validates Rust v1 implementation against Python v0 reference outputs.
  // Inputs (gm_mugration_inputs.json) define dataset/attribute combinations and
  // the exact capture parameters shared by the Python oracle and Rust replay.
  // Outputs (gm_mugration_outputs.json) were captured from v0 using gm_mugration_capture.
  //
  // Both v0 and v1 use iterative GTR inference (5 iterations of rate matrix
  // re-estimation followed by rate optimization via Brent's method).

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country(          "zika_20_country")]
  #[case::zika_20_country_weights(  "zika_20_country_weights")]
  #[trace]
  fn test_gm_mugration_outputs(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs[case];
    let expected = &outputs[case];
    let actual = run_gm_mugration_case(input)?;

    let expected_states = expected.states.clone();
    let actual_states = actual.gtr.states.clone();
    assert_eq!(expected_states, actual_states);

    let expected_n_states = expected.states.len();
    let actual_n_states = actual.gtr.n_states;
    assert_eq!(expected_n_states, actual_n_states);

    let expected_trait_assignments = expected.trait_assignments.clone();
    let actual_trait_assignments: BTreeMap<String, String> =
      actual.trait_assignments().iter().map(|(k, v)| (k.clone(), v.clone())).collect();
    assert_eq!(expected_trait_assignments, actual_trait_assignments);

    Ok(())
  }

  // Remaining 5 datasets diverge from v0 at ambiguous internal nodes.
  // Causes: D1 (apply_pseudo_counts on initial pi), D2 (root state uniform-threshold
  // filtering). Both are intentional v1 improvements over v0.
  #[rustfmt::skip]
  #[rstest]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
  #[case::dengue_20_country(        "dengue_20_country")]
  #[case::tb_20_cluster(            "tb_20_cluster")]
  #[case::rsv_a_20_country(         "rsv_a_20_country")]
  #[case::mpox_clade_ii_20_country( "mpox_clade_ii_20_country")]
  #[trace]
  #[ignore = "v0 parity: intentional v1 improvements (pseudo-count pi, root-state filtering)"]
  fn test_gm_mugration_outputs_v1_divergence(#[case] case: &str) -> Result<(), Report> {
    test_gm_mugration_outputs(case)
  }

  // Confidence profile comparison against v0 oracle.
  // Max observed error: 1.34e-2 at NODE_0000008 due to intentional v1
  // improvements (D1: pseudo-count pi, D2: root-state filtering).
  // See kb/issues/M-mugration-iterative-gtr.md.
  #[test]
  #[ignore = "v0 parity: max 1.34e-2 confidence divergence at ambiguous node (kb/issues/M-mugration-iterative-gtr.md)"]
  fn test_gm_mugration_confidence_zika() -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs["zika_20_country"];
    let expected = &outputs["zika_20_country"];
    let actual = run_gm_mugration_case(input)?;

    assert_eq!(expected.states, actual.confidence.states);

    for (node_name, expected_profile) in &expected.confidence {
      let actual_profile = actual
        .confidence
        .rows
        .iter()
        .find(|row| row.node == *node_name)
        .unwrap_or_else(|| panic!("missing confidence for node '{node_name}'"));
      let expected_arr = Array1::from_vec(expected_profile.clone());
      assert_abs_diff_eq!(expected_arr, actual_profile.profile, epsilon = 1e-6);
    }

    Ok(())
  }

  // Remaining confidence tests: max errors exceed 1e-2 at multiple nodes
  // due to D1/D2 intentional improvements.
  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country_weights(  "zika_20_country_weights")]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
  #[case::dengue_20_country(        "dengue_20_country")]
  #[case::tb_20_cluster(            "tb_20_cluster")]
  #[case::rsv_a_20_country(         "rsv_a_20_country")]
  #[case::mpox_clade_ii_20_country( "mpox_clade_ii_20_country")]
  #[trace]
  #[ignore = "v0 parity: confidence profiles exceed 1e-2 at multiple nodes (D1/D2 divergence)"]
  fn test_gm_mugration_confidence_outputs(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs[case];
    let expected = &outputs[case];
    let actual = run_gm_mugration_case(input)?;

    assert_eq!(expected.states, actual.confidence.states);

    for (node_name, expected_profile) in &expected.confidence {
      let actual_profile = actual
        .confidence
        .rows
        .iter()
        .find(|row| row.node == *node_name)
        .unwrap_or_else(|| panic!("missing confidence for node '{node_name}'"));
      let expected_arr = Array1::from_vec(expected_profile.clone());
      assert_abs_diff_eq!(expected_arr, actual_profile.profile, epsilon = 1e-10);
    }

    Ok(())
  }

  mod helpers {
    use crate::mugration::input::MugrationInput;
    use crate::mugration::mugration::execute_mugration;
    use crate::mugration::result::MugrationResult;
    use eyre::Report;
    use indexmap::IndexMap;
    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::path::PathBuf;
    use treetime_io::discrete_states_csv::read_discrete_attrs;
    use treetime_io::nwk::nwk_read_file;
    use treetime_utils::io::json::json_read_file;

    #[derive(Debug, Deserialize)]
    pub struct GmMugrationInput {
      pub tree_path: String,
      pub metadata_path: String,
      pub attribute: String,
      pub name_column: Option<String>,
      pub parameters: GmMugrationParameters,
    }

    #[derive(Debug, Deserialize)]
    pub struct GmMugrationParameters {
      pub missing_data: String,
      pub pc: Option<f64>,
      pub sampling_bias_correction: Option<f64>,
      pub weights_path: Option<String>,
      #[allow(dead_code)]
      pub verbose: usize,
      pub iterations: usize,
      #[allow(dead_code)]
      pub rng_seed: u64,
    }

    #[derive(Debug, Deserialize)]
    #[allow(dead_code)]
    pub struct GmMugrationOutput {
      pub states: Vec<String>,
      pub trait_assignments: BTreeMap<String, String>,
      pub confidence: BTreeMap<String, Vec<f64>>,
    }

    pub fn load_gm_mugration_inputs() -> IndexMap<String, GmMugrationInput> {
      let path = format!(
        "{}/src/commands/mugration/__tests__/__fixtures__/gm_mugration_inputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      json_read_file(&path).unwrap()
    }

    pub fn load_gm_mugration_outputs() -> IndexMap<String, GmMugrationOutput> {
      let path = format!(
        "{}/src/commands/mugration/__tests__/__fixtures__/gm_mugration_outputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      json_read_file(&path).unwrap()
    }

    pub fn run_gm_mugration_case(fixture: &GmMugrationInput) -> Result<MugrationResult, Report> {
      let project_root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");

      // Read tree directly
      let tree_path = project_root.join(&fixture.tree_path);
      let graph = nwk_read_file(&tree_path)?;

      // Read trait values using in-memory parsing
      let metadata_path = project_root.join(&fixture.metadata_path);
      let (attr_values, _attr_name) = read_discrete_attrs::<String>(
        &metadata_path,
        &fixture.name_column,
        &Some(fixture.attribute.clone()),
        |s| Ok(s.to_owned()),
      )?;
      let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

      // Read weights if provided
      let weights = match &fixture.parameters.weights_path {
        Some(weights_path) => {
          let weights_filepath = project_root.join(weights_path);
          let (map, _) = read_discrete_attrs::<f64>(
            &weights_filepath,
            &Some(fixture.attribute.clone()),
            &Some("weight".to_owned()),
            |s| Ok(s.parse::<f64>()?),
          )?;
          Some(map.into_iter().collect())
        },
        None => None,
      };

      // Build MugrationInput and execute
      let input = MugrationInput {
        graph,
        traits,
        attribute: fixture.attribute.clone(),
        weights,
        missing_data: fixture.parameters.missing_data.clone(),
        pc: fixture.parameters.pc,
        missing_weights_threshold: 0.5,
        iterations: fixture.parameters.iterations,
        sampling_bias_correction: fixture.parameters.sampling_bias_correction,
      };

      execute_mugration(input)
    }
  }
}
