#[cfg(test)]
mod tests {
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  use helpers::{load_gm_mugration_inputs, load_gm_mugration_outputs, run_gm_mugration_case};

  // Golden master tests for mugration discrete trait reconstruction.
  //
  // Validates Rust v1 implementation against Python v0 reference outputs.
  // Inputs (gm_mugration_inputs.json) define dataset/attribute combinations and
  // the exact capture parameters shared by the Python oracle and Rust replay.
  // Outputs (gm_mugration_outputs.json) were captured from v0 using gm_mugration_capture.
  //
  // v0 uses iterative GTR inference (5 iterations of rate matrix re-estimation) while
  // v1 uses uniform rates without iterative GTR. For datasets with strong phylogeographic
  // signal, the reconstructed states agree. Confidence profiles still diverge because the
  // one-shot v1 model remains sharper than the iteratively refined v0 model.
  //
  // Cases that pass: datasets where the phylogeographic signal overwhelms the GTR
  // model differences (zika, lassa).
  //
  // Cases that are ignored: datasets where iterative GTR re-estimation changes
  // the argmax at multiple internal nodes (dengue, tb, rsv, mpox). These will
  // pass once v1 implements iterative GTR optimization.
  //
  // Floating-point coverage:
  // - trait_assignments: discrete labels, exact match required
  // - confidence profiles: full-precision floats in oracle, compared with tight tolerance
  //   Blocked by missing iterative GTR; explicit ignored test preserves oracle contract.
  // - GTR equilibrium frequencies: not captured in oracle; v1 uses uniform frequencies
  //   while v0 iteratively refines from data. Adding coverage requires iterative GTR first.
  // - GTR rate parameters: not captured; same iterative GTR dependency.

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country(          "zika_20_country")]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
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
    let actual_trait_assignments: std::collections::BTreeMap<String, String> =
      actual.trait_assignments().iter().map(|(k, v)| (k.clone(), v.clone())).collect();
    assert_eq!(expected_trait_assignments, actual_trait_assignments);

    Ok(())
  }

  // Confidence profile comparison using full-precision floats.
  // The oracle stores floats at full precision; the test parses the CLI CSV output
  // back to floats for scientific comparison. This avoids double-rounding through
  // presentation formatting.
  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country(          "zika_20_country")]
  #[case::zika_20_country_weights(  "zika_20_country_weights")]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
  #[trace]
  #[ignore = "confidence parity depends on iterative GTR refinement: v0 re-estimates rates"]
  fn test_gm_mugration_confidence_outputs(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs[case];
    let expected = &outputs[case];
    let actual = run_gm_mugration_case(input)?;

    let expected_confidence_states = expected.states.clone();
    let actual_confidence_states = actual.confidence.states.clone();
    assert_eq!(expected_confidence_states, actual_confidence_states);

    for (node_name, expected_profile) in &expected.confidence {
      let actual_profile = actual
        .confidence
        .rows
        .iter()
        .find(|row| row.node == *node_name)
        .unwrap_or_else(|| panic!("missing confidence for node '{node_name}'"));
      assert_eq!(
        expected_profile.len(),
        actual_profile.profile.len(),
        "profile length mismatch for node '{node_name}'"
      );
      for (expected_val, actual_val) in expected_profile.iter().zip(actual_profile.profile.iter()) {
        assert_abs_diff_eq!(expected_val, actual_val, epsilon = 1e-10);
      }
    }

    Ok(())
  }

  // Ignored: v1 lacks iterative GTR optimization. v0 re-estimates the rate matrix
  // even when equilibrium frequencies are fixed from `--weights`, which still shifts
  // the argmax at ambiguous internal nodes. These tests will pass once v1 implements
  // iterative GTR. The captured v0 outputs remain valid oracles.
  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country_weights(  "zika_20_country_weights")]
  #[case::dengue_20_country(        "dengue_20_country")]
  #[case::tb_20_cluster(            "tb_20_cluster")]
  #[case::rsv_a_20_country(         "rsv_a_20_country")]
  #[case::mpox_clade_ii_20_country( "mpox_clade_ii_20_country")]
  #[trace]
  #[ignore = "v1 lacks iterative GTR optimization: v0 re-estimates equilibrium frequencies from data"]
  fn test_gm_mugration_outputs_iterative_gtr(#[case] case: &str) -> Result<(), Report> {
    test_gm_mugration_outputs(case)
  }

  mod helpers {
    use crate::commands::mugration::input::MugrationInput;
    use crate::commands::mugration::output::MugrationResult;
    use crate::commands::mugration::run::execute_mugration;
    use eyre::Report;
    use indexmap::IndexMap;
    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::path::PathBuf;
    use treetime_io::discrete_states_csv::read_discrete_attrs;
    use treetime_io::json::json_read_file;
    use treetime_io::nwk::nwk_read_file;

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
      #[allow(dead_code)]
      pub sampling_bias_correction: Option<f64>,
      pub weights_path: Option<String>,
      #[allow(dead_code)]
      pub verbose: usize,
      #[allow(dead_code)]
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
      };

      execute_mugration(input)
    }
  }
}
