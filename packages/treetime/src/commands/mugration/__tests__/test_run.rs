#[cfg(test)]
mod tests {
  use crate::commands::mugration::run::{
    apply_pseudo_counts, compute_pi_from_weights, compute_pi_uniform, run_mugration, validate_weight_coverage,
  };
  use crate::representation::discrete_states::DiscreteStates;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indexmap::IndexSet;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use std::fs;
  use treetime_utils::o;

  // ==========================================================================
  // Tests for extracted pure helpers (T11b-T11d)
  // ==========================================================================

  #[test]
  fn test_run_validate_weight_coverage_rejects_above_threshold() {
    let unique_values: IndexSet<String> = [o!("usa"), o!("germany"), o!("france"), o!("italy")]
      .into_iter()
      .collect();
    let weights_keys: IndexSet<String> = std::iter::once(o!("usa")).collect();
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
    let weights_keys: IndexSet<String> = std::iter::once(o!("usa")).collect();
    let missing_data = "?";
    let threshold = 0.5;

    let coverage = validate_weight_coverage(&unique_values, &weights_keys, missing_data, threshold).unwrap();
    let expected_missing: IndexSet<String> = std::iter::once(o!("germany")).collect();
    assert_eq!(expected_missing, coverage.missing_values);
    assert_abs_diff_eq!(0.5, coverage.missing_ratio, epsilon = 1e-10);
  }

  #[test]
  fn test_run_validate_weight_coverage_excludes_missing_data_marker() {
    let unique_values: IndexSet<String> = [o!("usa"), o!("?")].into_iter().collect();
    let weights_keys: IndexSet<String> = std::iter::once(o!("usa")).collect();
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

  // ==========================================================================
  // End-to-end file-based smoke tests (existing)
  // ==========================================================================

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
