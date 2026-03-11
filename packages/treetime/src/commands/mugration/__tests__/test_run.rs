#[cfg(test)]
mod tests {
  use crate::commands::mugration::run::run_mugration;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::fs;

  #[test]
  fn test_run_writes_expected_outputs_with_confidence() -> Result<(), Report> {
    let test_case = helpers::create_test_case("with-confidence")?;
    let args = helpers::make_args(&test_case, true);

    run_mugration(&args)?;

    helpers::assert_non_empty_file(&test_case.annotated_tree_path)?;
    helpers::assert_non_empty_file(&test_case.gtr_path)?;
    helpers::assert_non_empty_file(&test_case.traits_path)?;
    helpers::assert_non_empty_file(&test_case.confidence_path)?;

    let annotated_tree = fs::read_to_string(&test_case.annotated_tree_path)?;
    let expected_annotated_tree = helpers::expected_annotated_tree();
    assert_eq!(expected_annotated_tree, annotated_tree);

    let gtr = helpers::read_gtr_output(&test_case.gtr_path)?;
    let expected_gtr = helpers::expected_gtr_output();
    assert_eq!(expected_gtr.attribute, gtr.attribute);
    assert_eq!(expected_gtr.n_states, gtr.n_states);
    assert_eq!(expected_gtr.states, gtr.states);
    assert_eq!(expected_gtr.pi.len(), gtr.pi.len());
    assert_abs_diff_eq!(expected_gtr.pi[0], gtr.pi[0], epsilon = 1e-12);
    assert_abs_diff_eq!(expected_gtr.pi[1], gtr.pi[1], epsilon = 1e-12);
    assert_abs_diff_eq!(expected_gtr.mu, gtr.mu, epsilon = 1e-12);

    let trait_csv = fs::read_to_string(&test_case.traits_path)?;
    let expected_trait_csv = helpers::expected_trait_csv();
    assert_eq!(expected_trait_csv, trait_csv);

    let confidence_csv = fs::read_to_string(&test_case.confidence_path)?;
    let expected_confidence_csv = helpers::expected_confidence_csv();
    assert_eq!(expected_confidence_csv, confidence_csv);

    Ok(())
  }

  #[test]
  fn test_run_skips_confidence_output_when_not_requested() -> Result<(), Report> {
    let test_case = helpers::create_test_case("without-confidence")?;
    let args = helpers::make_args(&test_case, false);

    run_mugration(&args)?;

    let annotated_tree = fs::read_to_string(&test_case.annotated_tree_path)?;
    let expected_annotated_tree = helpers::expected_annotated_tree();
    assert_eq!(expected_annotated_tree, annotated_tree);

    let gtr = helpers::read_gtr_output(&test_case.gtr_path)?;
    let expected_gtr = helpers::expected_gtr_output();
    assert_eq!(expected_gtr.attribute, gtr.attribute);
    assert_eq!(expected_gtr.n_states, gtr.n_states);
    assert_eq!(expected_gtr.states, gtr.states);
    assert_eq!(expected_gtr.pi.len(), gtr.pi.len());
    assert_abs_diff_eq!(expected_gtr.pi[0], gtr.pi[0], epsilon = 1e-12);
    assert_abs_diff_eq!(expected_gtr.pi[1], gtr.pi[1], epsilon = 1e-12);
    assert_abs_diff_eq!(expected_gtr.mu, gtr.mu, epsilon = 1e-12);

    let trait_csv = fs::read_to_string(&test_case.traits_path)?;
    let expected_trait_csv = helpers::expected_trait_csv();
    assert_eq!(expected_trait_csv, trait_csv);
    assert!(!test_case.confidence_path.exists());

    Ok(())
  }

  mod helpers {
    use crate::commands::mugration::args::TreetimeMugrationArgs;
    use crate::{o, vec_of_owned};
    use eyre::Report;
    use indoc::indoc;
    use itertools::Itertools;
    use serde::Deserialize;
    use std::fs;
    use std::path::{Path, PathBuf};
    use treetime_utils::make_report;

    pub(super) fn create_test_case(test_name: &str) -> Result<TestCasePaths, Report> {
      let workspace_tmp_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..").join("tmp");
      let base_dir = workspace_tmp_dir.join(format!("test-run-mugration-{test_name}-{}", std::process::id()));
      if base_dir.exists() {
        fs::remove_dir_all(&base_dir)?;
      }
      let input_dir = base_dir.join("input");
      let output_dir = base_dir.join("output");
      fs::create_dir_all(&input_dir)?;
      fs::create_dir_all(&output_dir)?;

      let tree_path = input_dir.join("tree.nwk");
      let states_path = input_dir.join("states.tsv");
      fs::write(&tree_path, "(A:0.1,B:0.2)root;\n")?;
      fs::write(&states_path, "#name\tcountry\nA\tusa\nB\tgermany\n")?;

      Ok(TestCasePaths {
        tree_path,
        states_path,
        annotated_tree_path: output_dir.join("annotated_tree.nexus"),
        gtr_path: output_dir.join("gtr.json"),
        traits_path: output_dir.join("traits.csv"),
        confidence_path: output_dir.join("confidence.csv"),
        output_dir,
      })
    }

    pub(super) fn make_args(test_case: &TestCasePaths, with_confidence: bool) -> TreetimeMugrationArgs {
      TreetimeMugrationArgs {
        tree: Some(test_case.tree_path.clone()),
        attribute: "country".to_owned(),
        states: test_case.states_path.clone(),
        weights: None,
        name_column: None,
        confidence: with_confidence.then_some(test_case.confidence_path.clone()),
        pc: None,
        missing_data: "?".to_owned(),
        missing_weights_threshold: 0.5,
        sampling_bias_correction: None,
        outdir: test_case.output_dir.clone(),
        seed: None,
      }
    }

    pub(super) fn assert_non_empty_file(path: &Path) -> Result<(), Report> {
      if !path.exists() {
        return Err(make_report!("Expected file '{}' to exist", path.display()));
      }
      let contents = fs::read(path)?;
      if contents.is_empty() {
        return Err(make_report!("Expected file '{}' to be non-empty", path.display()));
      }
      Ok(())
    }

    pub(super) fn read_gtr_output(path: &Path) -> Result<GtrOutput, Report> {
      let contents = fs::read_to_string(path)?;
      Ok(serde_json::from_str(&contents)?)
    }

    pub(super) fn expected_annotated_tree() -> String {
      o!(indoc! {r#"
        #NEXUS
        Begin Taxa;
          Dimensions NTax=2;
          TaxLabels A B;
        End;
        Begin Trees;
          Tree tree1=(A:0.1[&country="usa"],B:0.2[&country="germany"])root[&country="usa"];;
        End;
        
      "#})
    }

    pub(super) fn expected_trait_csv() -> String {
      o!(indoc! {r#"
        node,country
        root,usa
        A,usa
        B,germany
      "#})
    }

    pub(super) fn expected_confidence_csv() -> String {
      o!(indoc! {r#"
        node,germany,usa
        root,0.333888,0.666112
        A,0.000000,1.000000
        B,1.000000,0.000000
      "#})
    }

    pub(super) fn expected_gtr_output() -> GtrOutput {
      GtrOutput {
        attribute: o!("country"),
        n_states: 2,
        states: vec_of_owned!["germany", "usa"],
        pi: vec![0.5, 0.5],
        mu: 0.5,
      }
    }

    #[derive(Debug)]
    pub(super) struct TestCasePaths {
      pub tree_path: PathBuf,
      pub states_path: PathBuf,
      pub annotated_tree_path: PathBuf,
      pub gtr_path: PathBuf,
      pub traits_path: PathBuf,
      pub confidence_path: PathBuf,
      pub output_dir: PathBuf,
    }

    #[derive(Debug, Deserialize, PartialEq)]
    pub(super) struct GtrOutput {
      pub attribute: String,
      pub n_states: usize,
      pub states: Vec<String>,
      pub pi: Vec<f64>,
      pub mu: f64,
    }
  }
}
