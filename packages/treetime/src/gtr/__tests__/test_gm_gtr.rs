#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{
    F81Params, HKY85Params, JC69Params, K80Params, T92Params, TN93Params, f81, hky85, jc69, k80, t92, tn93,
  };
  use crate::gtr::gtr::GTR;
  use eyre::Report;
  use rstest::rstest;

  use helpers::{compare_gtr, load_gm_gtr_inputs, load_gm_gtr_outputs};

  #[rstest]
  #[case::default("default")]
  #[case::mu_0_5("mu_0.5")]
  #[case::mu_2_0("mu_2.0")]
  #[trace]
  fn test_gm_gtr_jc69(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_gtr_inputs();
    let outputs = load_gm_gtr_outputs();
    let input = &inputs.jc69[case];
    let expected = &outputs.jc69[case];
    let gtr = jc69(JC69Params {
      mu: input.mu,
      alphabet: AlphabetName::Nuc,
    })?;
    compare_gtr(&gtr, expected);
    Ok(())
  }

  #[rstest]
  #[case::default("default")]
  #[case::kappa_0_5("kappa_0.5")]
  #[case::kappa_1_0("kappa_1.0")]
  #[case::kappa_2_0("kappa_2.0")]
  #[case::kappa_5_0("kappa_5.0")]
  #[case::mu_0_5_kappa_2_0("mu_0.5_kappa_2.0")]
  #[trace]
  fn test_gm_gtr_k80(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_gtr_inputs();
    let outputs = load_gm_gtr_outputs();
    let input = &inputs.k80[case];
    let expected = &outputs.k80[case];
    let gtr = k80(K80Params {
      mu: input.mu,
      kappa: input.kappa,
      alphabet: AlphabetName::Nuc,
    })?;
    compare_gtr(&gtr, expected);
    Ok(())
  }

  #[rstest]
  #[case::default("default")]
  #[case::custom_pi("custom_pi")]
  #[case::asymmetric_pi("asymmetric_pi")]
  #[case::mu_2_0("mu_2.0")]
  #[trace]
  fn test_gm_gtr_f81(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_gtr_inputs();
    let outputs = load_gm_gtr_outputs();
    let input = &inputs.f81[case];
    let expected = &outputs.f81[case];
    let gtr = f81(F81Params {
      mu: input.mu,
      pi: input.pi.clone(),
      alphabet: AlphabetName::Nuc,
    })?;
    compare_gtr(&gtr, expected);
    Ok(())
  }

  #[rstest]
  #[case::default("default")]
  #[case::kappa_2_0("kappa_2.0")]
  #[case::custom_pi("custom_pi")]
  #[case::custom_all("custom_all")]
  #[trace]
  fn test_gm_gtr_hky85(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_gtr_inputs();
    let outputs = load_gm_gtr_outputs();
    let input = &inputs.hky85[case];
    let expected = &outputs.hky85[case];
    let gtr = hky85(HKY85Params {
      mu: input.mu,
      kappa: input.kappa,
      pi: input.pi.clone(),
      alphabet: AlphabetName::Nuc,
    })?;
    compare_gtr(&gtr, expected);
    Ok(())
  }

  #[rstest]
  #[case::default("default")]
  #[case::high_gc("high_gc")]
  #[case::low_gc("low_gc")]
  #[case::kappa_2_0("kappa_2.0")]
  #[case::custom_all("custom_all")]
  #[trace]
  fn test_gm_gtr_t92(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_gtr_inputs();
    let outputs = load_gm_gtr_outputs();
    let input = &inputs.t92[case];
    let expected = &outputs.t92[case];
    let gtr = t92(T92Params {
      mu: input.mu,
      kappa: input.kappa,
      pi_GC: input.pi_gc,
      alphabet: AlphabetName::Nuc,
    })?;
    compare_gtr(&gtr, expected);
    Ok(())
  }

  #[rstest]
  #[case::default("default")]
  #[case::kappa1_0_5("kappa1_0.5")]
  #[case::kappa2_2_0("kappa2_2.0")]
  #[case::both_kappa("both_kappa")]
  #[case::custom_pi("custom_pi")]
  #[case::full_custom("full_custom")]
  #[trace]
  fn test_gm_gtr_tn93(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_gtr_inputs();
    let outputs = load_gm_gtr_outputs();
    let input = &inputs.tn93[case];
    let expected = &outputs.tn93[case];
    let gtr = tn93(TN93Params {
      mu: input.mu,
      kappa1: input.kappa1,
      kappa2: input.kappa2,
      pi: input.pi.clone(),
      alphabet: AlphabetName::Nuc,
    })?;
    compare_gtr(&gtr, expected);
    Ok(())
  }

  #[rstest]
  #[case::uniform_rates("uniform_rates")]
  #[case::ts_tv_bias("ts_tv_bias")]
  #[case::asymmetric_pi("asymmetric_pi")]
  #[case::mu_scaling("mu_scaling")]
  #[case::random("random")]
  #[trace]
  fn test_gm_gtr_custom(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_gtr_inputs();
    let outputs = load_gm_gtr_outputs();
    let input = &inputs.custom[case];
    let expected = &outputs.custom[case];
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let n_states = alphabet.n_canonical();
    let gtr = GTR::builder()
      .n_states(n_states)
      .mu(input.mu)
      .W(input.w.clone())
      .pi(input.pi.clone())
      .build()?;
    compare_gtr(&gtr, expected);
    Ok(())
  }

  mod helpers {
    use crate::gtr::gtr::GTR;
    use approx::assert_abs_diff_eq;
    use indexmap::IndexMap;
    use ndarray::{Array1, Array2};
    use serde::Deserialize;
    use std::fs::read_to_string;
    use treetime_utils::array::ndarray::sorted;
    use treetime_utils::array::serde::{array1_from_vec, array2_from_vec, option_array1_from_vec};

    #[derive(Debug, Deserialize)]
    pub(super) struct JC69Input {
      pub(crate) mu: f64,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct K80Input {
      pub(crate) mu: f64,
      pub(crate) kappa: f64,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct F81Input {
      pub(crate) mu: f64,
      #[serde(default, deserialize_with = "option_array1_from_vec")]
      pub(crate) pi: Option<Array1<f64>>,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct HKY85Input {
      pub(crate) mu: f64,
      pub(crate) kappa: f64,
      #[serde(default, deserialize_with = "option_array1_from_vec")]
      pub(crate) pi: Option<Array1<f64>>,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct T92Input {
      pub(crate) mu: f64,
      pub(crate) kappa: f64,
      pub(crate) pi_gc: f64,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct TN93Input {
      pub(crate) mu: f64,
      pub(crate) kappa1: f64,
      pub(crate) kappa2: f64,
      #[serde(default, deserialize_with = "option_array1_from_vec")]
      pub(crate) pi: Option<Array1<f64>>,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct CustomInput {
      pub(crate) mu: f64,
      #[serde(deserialize_with = "array1_from_vec")]
      pub(crate) pi: Array1<f64>,
      #[serde(rename = "W", deserialize_with = "array2_from_vec")]
      pub(crate) w: Array2<f64>,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct AllInputs {
      pub(crate) jc69: IndexMap<String, JC69Input>,
      pub(crate) k80: IndexMap<String, K80Input>,
      pub(crate) f81: IndexMap<String, F81Input>,
      pub(crate) hky85: IndexMap<String, HKY85Input>,
      pub(crate) t92: IndexMap<String, T92Input>,
      pub(crate) tn93: IndexMap<String, TN93Input>,
      pub(crate) custom: IndexMap<String, CustomInput>,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct ExpQtEntry {
      time: f64,
      #[serde(rename = "expQt", deserialize_with = "array2_from_vec")]
      exp_qt: Array2<f64>,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct GtrOutput {
      mu: f64,
      #[serde(deserialize_with = "array1_from_vec")]
      pi: Array1<f64>,
      #[serde(rename = "W", deserialize_with = "array2_from_vec")]
      w: Array2<f64>,
      #[serde(rename = "Q", deserialize_with = "array2_from_vec")]
      q: Array2<f64>,
      #[serde(deserialize_with = "array1_from_vec")]
      eigenvals: Array1<f64>,
      #[serde(rename = "expQts")]
      exp_qts: Vec<ExpQtEntry>,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct AllOutputs {
      pub(crate) jc69: IndexMap<String, GtrOutput>,
      pub(crate) k80: IndexMap<String, GtrOutput>,
      pub(crate) f81: IndexMap<String, GtrOutput>,
      pub(crate) hky85: IndexMap<String, GtrOutput>,
      pub(crate) t92: IndexMap<String, GtrOutput>,
      pub(crate) tn93: IndexMap<String, GtrOutput>,
      pub(crate) custom: IndexMap<String, GtrOutput>,
    }

    pub(super) fn load_gm_gtr_inputs() -> AllInputs {
      let path = format!(
        "{}/src/gtr/__tests__/__fixtures__/gm_gtr_inputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      let content = read_to_string(&path).unwrap_or_else(|_| panic!("Failed to read: {path}"));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {path}: {e}"))
    }

    pub(super) fn load_gm_gtr_outputs() -> AllOutputs {
      let path = format!(
        "{}/src/gtr/__tests__/__fixtures__/gm_gtr_outputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      let content = read_to_string(&path).unwrap_or_else(|_| panic!("Failed to read: {path}"));
      serde_json::from_str(&content).unwrap_or_else(|e| panic!("Failed to parse {path}: {e}"))
    }

    pub(super) fn compare_gtr(gtr: &GTR, expected: &GtrOutput) {
      assert_abs_diff_eq!(gtr.mu, expected.mu, epsilon = 1e-14);
      assert_abs_diff_eq!(gtr.pi, expected.pi, epsilon = 1e-14);
      assert_abs_diff_eq!(gtr.W, expected.w, epsilon = 1e-14);
      assert_abs_diff_eq!(gtr.Q(), expected.q, epsilon = 1e-14);

      assert_abs_diff_eq!(sorted(&gtr.eigvals), sorted(&expected.eigenvals), epsilon = 1e-14);

      for entry in &expected.exp_qts {
        let actual = gtr.expQt(entry.time);
        assert_abs_diff_eq!(actual, entry.exp_qt, epsilon = 1e-14);
      }
    }
  }
}
