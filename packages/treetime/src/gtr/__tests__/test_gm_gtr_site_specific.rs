#[cfg(test)]
mod tests {
  use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
  use crate::gtr::infer_gtr::site_specific::{InferGtrSiteSpecificOptions, infer_gtr_site_specific_impl};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::prelude::*;
  use rstest::rstest;

  use crate::gtr::__tests__::site_specific_support::{simulate_counts, value_to_array2, value_to_array3};
  use helpers::{load_gm_inputs, load_gm_outputs};

  // Golden master tests for site-specific GTR models.
  //
  // Validates Rust v1 GTRSiteSpecific against Python v0 GTR_site_specific outputs.
  // Inputs (gm_gtr_site_specific_inputs.json) define test parameters.
  // Outputs (gm_gtr_site_specific_outputs.json) captured from v0 using gm_gtr_site_specific_capture.
  //
  // See also:
  // - https://en.wikipedia.org/wiki/Characterization_test

  #[rstest]
  #[case::two_site_uniform_skewed("two_site_uniform_skewed")]
  #[case::three_site_heterogeneous("three_site_heterogeneous")]
  #[case::single_site_uniform("single_site_uniform")]
  fn test_gm_gtr_site_specific_expqt(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_inputs();
    let outputs = load_gm_outputs();
    let input = &inputs.site_specific[case];
    let expected = &outputs.site_specific[case];

    let gtr = helpers::build_from_input(input)?;

    // Compare eigenvalues (per-site, sorted per column for ordering independence)
    let expected_eigvals = value_to_array2(&expected.eigenvals);
    for a in 0..gtr.seq_len {
      let mut actual_col: Vec<f64> = gtr.eigvals.column(a).to_vec();
      let mut expected_col: Vec<f64> = expected_eigvals.column(a).to_vec();
      actual_col.sort_by_key(|x| ordered_float::OrderedFloat(*x));
      expected_col.sort_by_key(|x| ordered_float::OrderedFloat(*x));
      for (a_val, e_val) in actual_col.iter().zip(expected_col.iter()) {
        assert_abs_diff_eq!(a_val, e_val, epsilon = 1e-10);
      }
    }

    // Compare expQt at each captured time point
    for entry in &expected.exp_qts {
      let actual = gtr.expQt_raw(entry.time);
      let expected_qt = value_to_array3(&entry.exp_qt);
      assert_abs_diff_eq!(actual, expected_qt, epsilon = 1e-10);
    }

    Ok(())
  }

  #[rstest]
  #[case::two_site_uniform_skewed("two_site_uniform_skewed")]
  #[case::three_site_heterogeneous("three_site_heterogeneous")]
  fn test_gm_gtr_site_specific_propagate_evolve(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_inputs();
    let outputs = load_gm_outputs();
    let input = &inputs.site_specific[case];
    let expected = &outputs.site_specific[case];

    let gtr = helpers::build_from_input(input)?;

    for (profile_name, profile_output) in &expected.profiles {
      let profile_data = &inputs.profiles[profile_name.as_str()];
      let profile = value_to_array2(profile_data);

      let expected_prop = value_to_array2(&profile_output.propagate_profile);
      let expected_evol = value_to_array2(&profile_output.evolve);

      let actual_prop = gtr.propagate_profile(&profile, 0.5, false)?;
      let actual_evol = gtr.evolve(&profile, 0.5, false)?;

      assert_abs_diff_eq!(actual_prop, expected_prop, epsilon = 1e-10);
      assert_abs_diff_eq!(actual_evol, expected_evol, epsilon = 1e-10);
    }

    Ok(())
  }

  #[test]
  fn test_gm_gtr_site_specific_infer() -> Result<(), Report> {
    let inputs = load_gm_inputs();
    let outputs = load_gm_outputs();
    let infer_input = &inputs.infer["two_site"];
    let expected = &outputs.infer["two_site"];

    let n_states = infer_input["n_states"].as_u64().unwrap() as usize;
    let seq_len = infer_input["seq_len"].as_u64().unwrap() as usize;
    let pi = value_to_array2(&infer_input["pi"]);
    let mu: Vec<f64> = serde_json::from_value(infer_input["mu"].clone()).unwrap();
    let mu = Array1::from_vec(mu);
    let W = value_to_array2(&infer_input["W"]);
    let total_time = infer_input["total_time"].as_f64().unwrap();
    let pc = infer_input["pc"].as_f64().unwrap();

    let gtr_ref = GTRSiteSpecific::new(GTRSiteSpecificParams {
      n_states,
      seq_len,
      mu,
      W: Some(W),
      pi,
      approximate: false,
    })?;

    let counts = simulate_counts(&gtr_ref, total_time);

    let result = infer_gtr_site_specific_impl(
      &counts,
      &InferGtrSiteSpecificOptions {
        n_states,
        pc,
        max_iter: 30,
        dp: 1e-5,
        ..Default::default()
      },
    )?;

    let expected_pi = value_to_array2(&expected.pi);
    let expected_W = value_to_array2(&expected.W);

    // Compare inferred pi against v0 oracle.
    // Tolerance 1e-3: iterative convergence differs between numpy and ndarray due
    // to floating-point arithmetic order. With dp=1e-5 and total_time=5000, the
    // solvers converge to nearby but not identical fixed points.
    assert_abs_diff_eq!(result.pi, expected_pi, epsilon = 1e-3);

    // Compare W ratios against v0 oracle
    let v1_ref = result.W[[0, 1]];
    let v0_ref = expected_W[[0, 1]];
    for i in 0..n_states {
      for j in (i + 1)..n_states {
        let v1_ratio = result.W[[i, j]] / v1_ref;
        let v0_ratio = expected_W[[i, j]] / v0_ref;
        assert_abs_diff_eq!(v1_ratio, v0_ratio, epsilon = 1e-3);
      }
    }

    Ok(())
  }

  /// Approximate-mode golden master: compare v1 interpolation against v0 interpolation.
  #[rstest]
  #[case::two_site_uniform_skewed("two_site_uniform_skewed")]
  #[case::three_site_heterogeneous("three_site_heterogeneous")]
  #[case::single_site_uniform("single_site_uniform")]
  fn test_gm_gtr_site_specific_approximate(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_inputs();
    let outputs = load_gm_outputs();
    let input = &inputs.site_specific[case];
    let expected = &outputs.approximate[case];

    let gtr = helpers::build_from_input_approximate(input)?;

    for entry in &expected.exp_qts {
      let actual = gtr.expQt(entry.time)?;
      let expected_qt = value_to_array3(&entry.exp_qt);
      // v0 and v1 interpolation grids are identical (same 61-point non-uniform grid),
      // so interpolated outputs should match closely. Tolerance 1e-8 accounts for
      // floating-point differences in grid construction between numpy and ndarray.
      assert_abs_diff_eq!(actual, expected_qt, epsilon = 1e-8);
    }

    Ok(())
  }

  mod helpers {
    use super::*;
    use indexmap::IndexMap;
    use serde::Deserialize;
    use std::fs::read_to_string;

    pub fn parse_input(input: &serde_json::Value) -> (usize, usize, Array1<f64>, Option<Array2<f64>>, Array2<f64>) {
      let n_states = input["n_states"].as_u64().unwrap() as usize;
      let seq_len = input["seq_len"].as_u64().unwrap() as usize;
      let mu: Vec<f64> = serde_json::from_value(input["mu"].clone()).unwrap();
      let mu = Array1::from_vec(mu);
      let pi = value_to_array2(&input["pi"]);
      let W: Option<Array2<f64>> = if input["W"].is_null() {
        None
      } else {
        Some(value_to_array2(&input["W"]))
      };
      (n_states, seq_len, mu, W, pi)
    }

    pub fn build_from_input(input: &serde_json::Value) -> Result<GTRSiteSpecific, Report> {
      let (n_states, seq_len, mu, W, pi) = parse_input(input);
      GTRSiteSpecific::new(GTRSiteSpecificParams {
        n_states,
        seq_len,
        mu,
        W,
        pi,
        approximate: false,
      })
    }

    pub fn build_from_input_approximate(input: &serde_json::Value) -> Result<GTRSiteSpecific, Report> {
      let (n_states, seq_len, mu, W, pi) = parse_input(input);
      GTRSiteSpecific::new(GTRSiteSpecificParams {
        n_states,
        seq_len,
        mu,
        W,
        pi,
        approximate: true,
      })
    }

    #[derive(Debug, Deserialize)]
    pub struct ExpQtEntry {
      pub time: f64,
      #[serde(rename = "expQt")]
      pub exp_qt: serde_json::Value,
    }

    #[derive(Debug, Deserialize)]
    pub struct ProfileOutput {
      pub propagate_profile: serde_json::Value,
      pub evolve: serde_json::Value,
    }

    #[derive(Debug, Deserialize)]
    pub struct GtrSiteSpecificOutput {
      pub eigenvals: serde_json::Value,
      #[serde(rename = "expQts")]
      pub exp_qts: Vec<ExpQtEntry>,
      pub profiles: IndexMap<String, ProfileOutput>,
    }

    #[derive(Debug, Deserialize)]
    pub struct InferOutput {
      #[serde(rename = "W")]
      pub W: serde_json::Value,
      pub pi: serde_json::Value,
    }

    #[derive(Debug, Deserialize)]
    pub struct ApproxOutput {
      #[serde(rename = "expQts")]
      pub exp_qts: Vec<ExpQtEntry>,
    }

    #[derive(Debug, Deserialize)]
    pub struct AllOutputs {
      pub site_specific: IndexMap<String, GtrSiteSpecificOutput>,
      pub approximate: IndexMap<String, ApproxOutput>,
      pub infer: IndexMap<String, InferOutput>,
    }

    #[derive(Debug, Deserialize)]
    pub struct AllInputs {
      pub site_specific: IndexMap<String, serde_json::Value>,
      pub profiles: IndexMap<String, serde_json::Value>,
      pub infer: IndexMap<String, serde_json::Value>,
    }

    pub fn load_gm_inputs() -> AllInputs {
      let path = format!(
        "{}/src/gtr/__tests__/__fixtures__/gm_gtr_site_specific_inputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      let content = read_to_string(&path).unwrap_or_else(|_| panic!("Failed to read: {path}"));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {path}: {e}"))
    }

    pub fn load_gm_outputs() -> AllOutputs {
      let path = format!(
        "{}/src/gtr/__tests__/__fixtures__/gm_gtr_site_specific_outputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      let content = read_to_string(&path).unwrap_or_else(|_| panic!("Failed to read: {path}"));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {path}: {e}"))
    }
  }
}
