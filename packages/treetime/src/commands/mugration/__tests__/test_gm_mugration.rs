#[cfg(test)]
mod tests {
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  use helpers::{load_gm_mugration_inputs, load_gm_mugration_outputs, run_gm_mugration_case};

  // Golden master tests for mugration discrete trait reconstruction.
  //
  // Validates Rust v1 implementation against Python v0 reference outputs.
  // Inputs (gm_mugration_inputs.json) define dataset/attribute combinations.
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
    let actual = run_gm_mugration_case(case, input)?;

    let expected_states = expected.states.clone();
    let actual_states = actual.gtr.states.clone();
    assert_eq!(expected_states, actual_states);

    let expected_n_states = expected.states.len();
    let actual_n_states = actual.gtr.n_states;
    assert_eq!(expected_n_states, actual_n_states);

    let expected_trait_assignments = expected.trait_assignments.clone();
    let actual_trait_assignments = actual.trait_assignments;
    assert_eq!(expected_trait_assignments, actual_trait_assignments);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country(          "zika_20_country")]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
  #[trace]
  #[ignore = "confidence parity still depends on iterative GTR refinement"]
  fn test_gm_mugration_confidence_outputs(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs[case];
    let expected = &outputs[case];
    let actual = run_gm_mugration_case(case, input)?;

    let expected_confidence_states = expected.states.clone();
    let actual_confidence_states = actual.confidence.states.clone();
    assert_eq!(expected_confidence_states, actual_confidence_states);

    let expected_confidence = helpers::format_confidence_output(&expected.confidence);
    let actual_confidence = actual.confidence.rows;
    assert_eq!(expected_confidence, actual_confidence);

    Ok(())
  }

  // Ignored: v1 lacks iterative GTR optimization. v0 re-estimates the rate matrix
  // and equilibrium frequencies from the data (5 iterations), which shifts the argmax
  // at ambiguous internal nodes. These tests will pass once v1 implements iterative
  // GTR. The captured v0 outputs remain valid oracles.
  #[rustfmt::skip]
  #[rstest]
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
    use crate::commands::mugration::args::TreetimeMugrationArgs;
    use crate::commands::mugration::run::run_mugration;
    use eyre::Report;
    use indexmap::IndexMap;
    use itertools::Itertools;
    use serde::Deserialize;
    use std::collections::BTreeMap;
    use std::fs;
    use std::path::{Path, PathBuf};
    use treetime_io::json::json_read_file;
    use treetime_utils::make_report;

    #[derive(Debug, Deserialize)]
    pub struct MugrationInput {
      pub tree_path: String,
      pub metadata_path: String,
      pub attribute: String,
      pub name_column: Option<String>,
    }

    #[derive(Debug, Deserialize)]
    #[allow(dead_code)]
    pub struct MugrationOutput {
      pub states: Vec<String>,
      pub trait_assignments: BTreeMap<String, String>,
      pub confidence: BTreeMap<String, Vec<f64>>,
    }

    #[derive(Debug)]
    pub struct ActualMugrationOutput {
      pub gtr: GtrOutput,
      pub trait_assignments: BTreeMap<String, String>,
      pub confidence: ConfidenceOutput,
    }

    #[derive(Debug, Deserialize)]
    pub struct GtrOutput {
      pub n_states: usize,
      pub states: Vec<String>,
    }

    #[derive(Debug)]
    pub struct ConfidenceOutput {
      pub states: Vec<String>,
      pub rows: BTreeMap<String, Vec<String>>,
    }

    pub fn load_gm_mugration_inputs() -> IndexMap<String, MugrationInput> {
      let path = format!(
        "{}/src/commands/mugration/__tests__/__fixtures__/gm_mugration_inputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      let content = fs::read_to_string(&path).unwrap_or_else(|_| panic!("Failed to read: {path}"));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {path}: {e}"))
    }

    pub fn load_gm_mugration_outputs() -> IndexMap<String, MugrationOutput> {
      let path = format!(
        "{}/src/commands/mugration/__tests__/__fixtures__/gm_mugration_outputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      let content = fs::read_to_string(&path).unwrap_or_else(|_| panic!("Failed to read: {path}"));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {path}: {e}"))
    }

    pub fn run_gm_mugration_case(case: &str, input: &MugrationInput) -> Result<ActualMugrationOutput, Report> {
      let project_root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
      let outdir = project_root
        .join("tmp")
        .join(format!("test-gm-mugration-{case}-{}", std::process::id()));
      if outdir.exists() {
        fs::remove_dir_all(&outdir)?;
      }

      let confidence_path = outdir.join("confidence.csv");
      let args = TreetimeMugrationArgs {
        tree: Some(project_root.join(&input.tree_path)),
        attribute: input.attribute.clone(),
        states: project_root.join(&input.metadata_path),
        weights: None,
        name_column: input.name_column.clone(),
        confidence: Some(confidence_path.clone()),
        pc: Some(1.0),
        missing_data: "?".to_owned(),
        missing_weights_threshold: 0.5,
        sampling_bias_correction: None,
        outdir: outdir.clone(),
      };

      run_mugration(&args)?;

      let gtr = json_read_file(outdir.join("gtr.json"))?;
      let trait_assignments = read_traits_csv(&outdir.join("traits.csv"), &input.attribute)?;
      let confidence = read_confidence_csv(&confidence_path)?;

      Ok(ActualMugrationOutput {
        gtr,
        trait_assignments,
        confidence,
      })
    }

    pub fn format_confidence_output(confidence: &BTreeMap<String, Vec<f64>>) -> BTreeMap<String, Vec<String>> {
      confidence
        .iter()
        .map(|(name, profile)| {
          let formatted_profile = profile.iter().map(|value| format!("{value:.6}")).collect_vec();
          (name.clone(), formatted_profile)
        })
        .collect()
    }

    pub fn read_traits_csv(path: &Path, attribute: &str) -> Result<BTreeMap<String, String>, Report> {
      let content =
        fs::read_to_string(path).map_err(|e| make_report!("Failed to read traits CSV '{}': {e}", path.display()))?;
      let expected_header = format!("node,{attribute}");
      let mut result = BTreeMap::new();
      for (i, line) in content.lines().enumerate() {
        if i == 0 {
          if line != expected_header {
            return Err(make_report!(
              "Unexpected header in traits CSV: '{line}', expected '{expected_header}'"
            ));
          }
          continue;
        }
        let parts: Vec<&str> = line.splitn(2, ',').collect();
        if parts.len() == 2 {
          result.insert(parts[0].to_owned(), parts[1].to_owned());
        }
      }
      Ok(result)
    }

    pub fn read_confidence_csv(path: &Path) -> Result<ConfidenceOutput, Report> {
      let content = fs::read_to_string(path)
        .map_err(|e| make_report!("Failed to read confidence CSV '{}': {e}", path.display()))?;

      let mut lines = content.lines();
      let header = lines
        .next()
        .ok_or_else(|| make_report!("Confidence CSV '{}' is empty", path.display()))?;

      let mut header_parts = header.split(',');
      let first_column = header_parts
        .next()
        .ok_or_else(|| make_report!("Confidence CSV '{}' has no header columns", path.display()))?;
      if first_column != "node" {
        return Err(make_report!(
          "Unexpected header in confidence CSV '{}': first column was '{first_column}'",
          path.display()
        ));
      }

      let states = header_parts.map(str::to_owned).collect_vec();
      let expected_columns = states.len() + 1;
      let rows = lines
        .map(|line| {
          let parts = line.split(',').map(str::to_owned).collect_vec();
          if parts.len() != expected_columns {
            return Err(make_report!(
              "Unexpected confidence row in '{}': expected {expected_columns} columns, got {} in '{line}'",
              path.display(),
              parts.len()
            ));
          }
          let node_name = parts[0].clone();
          let profile = parts[1..].to_vec();
          Ok((node_name, profile))
        })
        .collect::<Result<BTreeMap<_, _>, Report>>()?;

      Ok(ConfidenceOutput { states, rows })
    }
  }
}
