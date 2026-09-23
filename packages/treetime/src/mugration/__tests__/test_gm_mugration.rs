#[cfg(test)]
mod tests {
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  use helpers::{
    confidence_by_name, load_gm_mugration_inputs, load_gm_mugration_outputs, run_gm_mugration_case, states_vec,
    trait_assignments_by_name,
  };

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country(          "zika_20_country")]
  #[case::zika_20_country_weights(  "zika_20_country_weights")]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
  #[trace]
  fn test_gm_mugration_outputs(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs[case];
    let expected = &outputs[case];
    let (output, names) = run_gm_mugration_case(input)?;

    let expected_states = expected.states.clone();
    let actual_states = states_vec(&output);
    assert_eq!(expected_states, actual_states);

    let expected_n_states = expected.states.len();
    let actual_n_states = output.n_states;
    assert_eq!(expected_n_states, actual_n_states);

    let expected_trait_assignments = expected.trait_assignments.clone();
    let actual_trait_assignments = trait_assignments_by_name(&output, &names);
    assert_eq!(expected_trait_assignments, actual_trait_assignments);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dengue_20_country(        "dengue_20_country")]
  #[case::tb_20_cluster(            "tb_20_cluster")]
  #[case::rsv_a_20_country(         "rsv_a_20_country")]
  #[case::mpox_clade_ii_20_country( "mpox_clade_ii_20_country")]
  #[trace]
  #[ignore = "v0 parity: residual ~1e-3 marginal divergence tips argmax at ambiguous nodes (kb/issues/M-mugration-iterative-gtr.md)"]
  fn test_gm_mugration_outputs_v1_divergence(#[case] case: &str) -> Result<(), Report> {
    test_gm_mugration_outputs(case)
  }

  #[test]
  #[ignore = "v0 parity: residual ~1e-3 marginal-confidence divergence (kb/issues/M-mugration-iterative-gtr.md)"]
  fn test_gm_mugration_confidence_zika() -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs["zika_20_country"];
    let expected = &outputs["zika_20_country"];
    let (output, names) = run_gm_mugration_case(input)?;

    assert_eq!(expected.states, states_vec(&output));

    let actual_confidence = confidence_by_name(&output, &names);
    for (node_name, expected_profile) in &expected.confidence {
      let actual_profile = actual_confidence
        .get(node_name)
        .unwrap_or_else(|| panic!("missing confidence for node '{node_name}'"));
      let expected_arr = Array1::from_vec(expected_profile.clone());
      assert_abs_diff_eq!(expected_arr, actual_profile, epsilon = 1e-6);
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country_weights(  "zika_20_country_weights")]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
  #[case::dengue_20_country(        "dengue_20_country")]
  #[case::tb_20_cluster(            "tb_20_cluster")]
  #[case::rsv_a_20_country(         "rsv_a_20_country")]
  #[case::mpox_clade_ii_20_country( "mpox_clade_ii_20_country")]
  #[trace]
  #[ignore = "v0 parity: residual marginal-confidence divergence from v0 (kb/issues/M-mugration-iterative-gtr.md)"]
  fn test_gm_mugration_confidence_outputs(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs[case];
    let expected = &outputs[case];
    let (output, names) = run_gm_mugration_case(input)?;

    assert_eq!(expected.states, states_vec(&output));

    let actual_confidence = confidence_by_name(&output, &names);
    for (node_name, expected_profile) in &expected.confidence {
      let actual_profile = actual_confidence
        .get(node_name)
        .unwrap_or_else(|| panic!("missing confidence for node '{node_name}'"));
      let expected_arr = Array1::from_vec(expected_profile.clone());
      assert_abs_diff_eq!(expected_arr, actual_profile, epsilon = 1e-10);
    }

    Ok(())
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::mugration::pipeline::{MugrationInput, MugrationOutput, MugrationParams, run};
    use eyre::Report;
    use indexmap::IndexMap;
    use ndarray::Array1;
    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::path::PathBuf;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::csv::default_name_candidates;
    use treetime_io::discrete_states_csv::read_discrete_attrs;
    use treetime_io::nwk::nwk_read_file;
    use treetime_utils::io::json::json_read_file;

    #[derive(Debug, Deserialize)]
    pub struct GmMugrationInput {
      tree_path: String,
      metadata_path: String,
      attribute: String,
      name_column: Option<String>,
      parameters: GmMugrationParameters,
    }

    #[derive(Debug, Deserialize)]
    pub struct GmMugrationParameters {
      pub missing_data: String,
      pub pc: Option<f64>,
      pub sampling_bias_correction: Option<f64>,
      pub weights_path: Option<String>,
      pub iterations: usize,
    }

    #[derive(Debug, Deserialize)]
    pub struct GmMugrationOutput {
      pub states: Vec<String>,
      pub trait_assignments: BTreeMap<String, String>,
      pub confidence: BTreeMap<String, Vec<f64>>,
    }

    pub fn load_gm_mugration_inputs() -> IndexMap<String, GmMugrationInput> {
      let path = format!(
        "{}/src/mugration/__tests__/__fixtures__/gm_mugration_inputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      json_read_file(&path).unwrap()
    }

    pub fn load_gm_mugration_outputs() -> IndexMap<String, GmMugrationOutput> {
      let path = format!(
        "{}/src/mugration/__tests__/__fixtures__/gm_mugration_outputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      json_read_file(&path).unwrap()
    }

    pub fn run_gm_mugration_case(
      fixture: &GmMugrationInput,
    ) -> Result<(MugrationOutput, BTreeMap<GraphNodeKey, Option<String>>), Report> {
      let project_root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");

      let tree_path = project_root.join(&fixture.tree_path);
      let nwk_parsed = nwk_read_file(&tree_path)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let metadata_path = project_root.join(&fixture.metadata_path);
      let (attr_values, _attr_name) = read_discrete_attrs::<String>(
        &metadata_path,
        &[',', '\t', ';'],
        &default_name_candidates(),
        &fixture.name_column,
        &Some(fixture.attribute.clone()),
        |s| Ok(s.to_owned()),
      )?;
      let traits: BTreeMap<String, String> = attr_values.into_iter().collect();

      let weights = match &fixture.parameters.weights_path {
        Some(weights_path) => {
          let weights_filepath = project_root.join(weights_path);
          let (map, _) = read_discrete_attrs::<f64>(
            &weights_filepath,
            &[',', '\t', ';'],
            &[],
            &Some(fixture.attribute.clone()),
            &Some("weight".to_owned()),
            |s| Ok(s.parse::<f64>()?),
          )?;
          Some(map.into_iter().collect())
        },
        None => None,
      };

      let params = MugrationParams {
        missing_data: fixture.parameters.missing_data.clone(),
        pc: fixture.parameters.pc,
        missing_weights_threshold: 0.5,
        iterations: fixture.parameters.iterations,
        sampling_bias_correction: fixture.parameters.sampling_bias_correction,
        smooth_initial_pi: false,
        filter_uninformative_root: false,
      };
      let input = MugrationInput {
        graph,
        traits,
        weights,
        branch_lengths,
      };
      let output = run(&params, input, &names, &NoopCancel).map_err(|err| err.into_report())?;
      Ok((output, names))
    }

    pub fn states_vec(output: &MugrationOutput) -> Vec<String> {
      output.states.iter().map(|s| s.to_owned()).collect()
    }

    pub fn trait_assignments_by_name(
      output: &MugrationOutput,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
    ) -> BTreeMap<String, String> {
      output
        .graph
        .get_nodes()
        .filter_map(|node| {
          let key = node.key();
          let name = names[&key].clone().unwrap_or_else(|| format!("node_{}", key.0));
          output.reconstructed_traits[&key].clone().map(|value| (name, value))
        })
        .collect()
    }

    pub fn confidence_by_name(
      output: &MugrationOutput,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
    ) -> BTreeMap<String, Array1<f64>> {
      output
        .graph
        .get_nodes()
        .filter_map(|node| {
          let key = node.key();
          let name = names[&key].clone().unwrap_or_else(|| format!("node_{}", key.0));
          output.confidences[&key].clone().map(|profile| (name, profile))
        })
        .collect()
    }
  }
}
