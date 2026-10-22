#[cfg(test)]
mod tests {
  use crate::commands::shared::argument_error::ArgumentError;
  use crate::commands::shared::model::{GtrModelNameCli, parse_model_params};
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::{Array1, Array2, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use strum::IntoEnumIterator;
  use treetime::alphabet::alphabet::AlphabetName;
  use treetime::gtr::get_gtr::{GtrModel, HKY85Params, Jtt92Params, K80Params, TN93Params};
  use treetime_utils::io::json::{JsonPretty, json_write_str};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::jc69_default(        GtrModelNameCli::JC69,  &[],                                                     ("jc69",  "default"))]
  #[case::jc69_mu(             GtrModelNameCli::JC69,  &["mu=0.5"],                                             ("jc69",  "mu_0.5"))]
  #[case::jc69_mu_2(           GtrModelNameCli::JC69,  &["mu=2.0"],                                             ("jc69",  "mu_2.0"))]
  #[case::k80_default(         GtrModelNameCli::K80,   &[],                                                     ("k80",   "default"))]
  #[case::k80_kappa_low(       GtrModelNameCli::K80,   &["kappa=0.5"],                                          ("k80",   "kappa_0.5"))]
  #[case::k80_kappa_one(       GtrModelNameCli::K80,   &["kappa=1.0"],                                          ("k80",   "kappa_1.0"))]
  #[case::k80_kappa_high(      GtrModelNameCli::K80,   &["kappa=5.0"],                                          ("k80",   "kappa_5.0"))]
  #[case::k80_mu_kappa(        GtrModelNameCli::K80,   &["mu=0.5", "kappa=2.0"],                                ("k80",   "mu_0.5_kappa_2.0"))]
  #[case::f81_default(         GtrModelNameCli::F81,   &[],                                                     ("f81",   "default"))]
  #[case::f81_pi(              GtrModelNameCli::F81,   &["pi=0.1,0.2,0.3,0.4"],                                 ("f81",   "custom_pi"))]
  #[case::f81_pis(             GtrModelNameCli::F81,   &["pis=0.5,0.1,0.1,0.3"],                                ("f81",   "asymmetric_pi"))]
  #[case::f81_mu(              GtrModelNameCli::F81,   &["mu=2.0"],                                             ("f81",   "mu_2.0"))]
  #[case::hky85_default(       GtrModelNameCli::HKY85, &[],                                                     ("hky85", "default"))]
  #[case::hky85_kappa(         GtrModelNameCli::HKY85, &["kappa=2.0"],                                          ("hky85", "kappa_2.0"))]
  #[case::hky85_capital_pi(    GtrModelNameCli::HKY85, &["Pi=0.1,0.2,0.3,0.4"],                                 ("hky85", "custom_pi"))]
  #[case::hky85_all(           GtrModelNameCli::HKY85, &["mu=1.5", "kappa=3.0", "Pis=0.3,0.2,0.2,0.3"],         ("hky85", "custom_all"))]
  #[case::t92_default(         GtrModelNameCli::T92,   &[],                                                     ("t92",   "default"))]
  #[case::t92_high_gc(         GtrModelNameCli::T92,   &["pi_GC=0.7"],                                          ("t92",   "high_gc"))]
  #[case::t92_low_gc(          GtrModelNameCli::T92,   &["pi_GC=0.3"],                                          ("t92",   "low_gc"))]
  #[case::t92_kappa(           GtrModelNameCli::T92,   &["kappa=2.0"],                                          ("t92",   "kappa_2.0"))]
  #[case::t92_all(             GtrModelNameCli::T92,   &["mu=2.0", "kappa=3.0", "pi_GC=0.65"],                  ("t92",   "custom_all"))]
  #[case::tn93_default(        GtrModelNameCli::TN93,  &[],                                                     ("tn93",  "default"))]
  #[case::tn93_kappa1(         GtrModelNameCli::TN93,  &["kappa1=0.5"],                                         ("tn93",  "kappa1_0.5"))]
  #[case::tn93_kappa2(         GtrModelNameCli::TN93,  &["kappa2=2.0"],                                         ("tn93",  "kappa2_2.0"))]
  #[case::tn93_both_kappa(     GtrModelNameCli::TN93,  &["kappa1=0.3", "kappa2=2.5"],                           ("tn93",  "both_kappa"))]
  #[case::tn93_pi(             GtrModelNameCli::TN93,  &["pi=0.3,0.2,0.25,0.25"],                               ("tn93",  "custom_pi"))]
  #[case::tn93_all(            GtrModelNameCli::TN93,  &["mu=1.5", "kappa1=0.5", "kappa2=3.0", "pi=0.2,0.3,0.3,0.2"], ("tn93", "full_custom"))]
  #[trace]
  fn test_model_params_build_the_v0_model(
    #[case] model: GtrModelNameCli,
    #[case] items: &[&str],
    #[case] (fixture_model, fixture_case): (&str, &str),
  ) -> Result<(), Report> {
    let expected = &helpers::v0_models()?[fixture_model][fixture_case];

    let gtr = parse_model_params(model, &helpers::owned(items))?.build(AlphabetName::Nuc)?;

    pretty_assert_abs_diff_eq!(expected.mu, gtr.mu(), epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(expected.pi.clone(), gtr.pi().clone(), epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(expected.w.clone(), gtr.W().clone(), epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::jc69_mu_zero(  GtrModelNameCli::JC69,  &["mu=0"],    (0.0, array![
    [0.0,       4.0 / 3.0, 4.0 / 3.0, 4.0 / 3.0],
    [4.0 / 3.0, 0.0,       4.0 / 3.0, 4.0 / 3.0],
    [4.0 / 3.0, 4.0 / 3.0, 0.0,       4.0 / 3.0],
    [4.0 / 3.0, 4.0 / 3.0, 4.0 / 3.0, 0.0      ],
  ]))]
  #[case::k80_kappa_zero(GtrModelNameCli::K80,   &["kappa=0"], (1.0, array![
    [0.0, 2.0, 0.0, 2.0],
    [2.0, 0.0, 2.0, 0.0],
    [0.0, 2.0, 0.0, 2.0],
    [2.0, 0.0, 2.0, 0.0],
  ]))]
  #[trace]
  fn test_model_params_zero_rates_build_the_model_definition(
    #[case] model: GtrModelNameCli,
    #[case] items: &[&str],
    #[case] (mu, w): (f64, Array2<f64>),
  ) -> Result<(), Report> {
    let gtr = parse_model_params(model, &helpers::owned(items))?.build(AlphabetName::Nuc)?;

    pretty_assert_abs_diff_eq!(mu, gtr.mu(), epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(array![0.25, 0.25, 0.25, 0.25], gtr.pi().clone(), epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(w, gtr.W().clone(), epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::infer(      GtrModelNameCli::Infer, &[],                                  GtrModel::Infer)]
  #[case::jtt92_mu(   GtrModelNameCli::Jtt92, &["mu=0.5"],                          GtrModel::Jtt92(Jtt92Params { mu: 0.5, ..Jtt92Params::default() }))]
  #[case::hky85_help( GtrModelNameCli::HKY85, &["kappa=0.2", "pis=0.25,0.25,0.25,0.25"], GtrModel::HKY85(HKY85Params { kappa: 0.2, pi: Some(array![0.25, 0.25, 0.25, 0.25]), ..HKY85Params::default() }))]
  #[case::kappa_zero( GtrModelNameCli::K80,   &["kappa=0"],                         GtrModel::K80(K80Params { kappa: 0.0, ..K80Params::default() }))]
  #[case::mu_zero(    GtrModelNameCli::TN93,  &["mu=0", "kappa1=0", "kappa2=0"],    GtrModel::TN93(TN93Params { mu: 0.0, kappa1: 0.0, kappa2: 0.0, ..TN93Params::default() }))]
  #[case::exponent(   GtrModelNameCli::HKY85, &["kappa=2e-1"],                      GtrModel::HKY85(HKY85Params { kappa: 0.2, ..HKY85Params::default() }))]
  #[trace]
  fn test_model_params_parse_into_the_core_value(
    #[case] model: GtrModelNameCli,
    #[case] items: &[&str],
    #[case] expected: GtrModel,
  ) -> Result<(), Report> {
    assert_eq!(expected, parse_model_params(model, &helpers::owned(items))?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::infer_with_params( GtrModelNameCli::Infer, &["kappa=0.2"],               ("--model-params needs a named --model: infer fits the model from the data and takes no parameters, but got 'kappa=0.2'", &["model", "model_params"][..]))]
  #[case::no_equals(         GtrModelNameCli::K80,   &["kappa"],                   ("--model-params takes key=value items, but got 'kappa'", &["model_params"][..]))]
  #[case::empty_key(         GtrModelNameCli::K80,   &["=0.2"],                    ("--model-params takes key=value items, but got '=0.2'", &["model_params"][..]))]
  #[case::empty_value(       GtrModelNameCli::K80,   &["kappa="],                  ("--model-params takes key=value items, but got 'kappa='", &["model_params"][..]))]
  #[case::two_equals(        GtrModelNameCli::K80,   &["kappa=0.2=0.3"],           ("--model-params takes key=value items, but got 'kappa=0.2=0.3'", &["model_params"][..]))]
  #[case::key_of_other_model(GtrModelNameCli::JC69,  &["kappa=0.2"],               ("--model-params item 'kappa=0.2': the model jc69 has no parameter 'kappa'; its parameters are mu", &["model", "model_params"][..]))]
  #[case::pi_on_k80(         GtrModelNameCli::K80,   &["pi=0.25,0.25,0.25,0.25"],  ("--model-params item 'pi=0.25,0.25,0.25,0.25': the model k80 has no parameter 'pi'; its parameters are mu, kappa", &["model", "model_params"][..]))]
  #[case::alphabet(          GtrModelNameCli::Jtt92, &["alphabet=aa"],             ("--model-params item 'alphabet=aa': the model jtt92 has no parameter 'alphabet'; its parameters are mu", &["model", "model_params"][..]))]
  #[case::unknown_key(       GtrModelNameCli::TN93,  &["kappa=2"],                 ("--model-params item 'kappa=2': the model tn93 has no parameter 'kappa'; its parameters are mu, kappa1, kappa2, pi", &["model", "model_params"][..]))]
  #[case::key_case(          GtrModelNameCli::T92,   &["pi_gc=0.4"],               ("--model-params item 'pi_gc=0.4': the model t92 has no parameter 'pi_gc'; its parameters are mu, kappa, pi_GC", &["model", "model_params"][..]))]
  #[case::repeated_key(      GtrModelNameCli::K80,   &["kappa=0.2", "kappa=0.3"],  ("--model-params sets kappa twice: 'kappa=0.2' and 'kappa=0.3'", &["model_params"][..]))]
  #[case::repeated_alias(    GtrModelNameCli::F81,   &["pi=0.1,0.2,0.3,0.4", "Pis=0.4,0.3,0.2,0.1"], ("--model-params sets pi twice: 'pi=0.1,0.2,0.3,0.4' and 'Pis=0.4,0.3,0.2,0.1'", &["model_params"][..]))]
  #[case::not_a_number(      GtrModelNameCli::K80,   &["kappa=high"],              ("--model-params item 'kappa=high': kappa must be a non-negative finite number", &["model_params"][..]))]
  #[case::negative_kappa1(   GtrModelNameCli::TN93,  &["kappa1=-0.5"],             ("--model-params item 'kappa1=-0.5': kappa1 must be a non-negative finite number", &["model_params"][..]))]
  #[case::negative_mu(       GtrModelNameCli::JC69,  &["mu=-1"],                   ("--model-params item 'mu=-1': mu must be a non-negative finite number", &["model_params"][..]))]
  #[case::infinite_kappa2(   GtrModelNameCli::TN93,  &["kappa2=inf"],              ("--model-params item 'kappa2=inf': kappa2 must be a non-negative finite number", &["model_params"][..]))]
  #[case::gc_at_one(         GtrModelNameCli::T92,   &["pi_GC=1"],                 ("--model-params item 'pi_GC=1': pi_GC must be a number between 0 and 1, exclusive", &["model_params"][..]))]
  #[case::gc_at_zero(        GtrModelNameCli::T92,   &["pi_GC=0"],                 ("--model-params item 'pi_GC=0': pi_GC must be a number between 0 and 1, exclusive", &["model_params"][..]))]
  #[case::pi_zero(           GtrModelNameCli::TN93,  &["pi=0,0.5,0.25,0.25"],      ("--model-params item 'pi=0,0.5,0.25,0.25': each frequency of pi must be a positive finite number", &["model_params"][..]))]
  #[case::pi_empty_entry(    GtrModelNameCli::TN93,  &["pi=0.5,,0.25,0.25"],       ("--model-params item 'pi=0.5,,0.25,0.25': each frequency of pi must be a positive finite number", &["model_params"][..]))]
  #[trace]
  fn test_model_params_reject_invalid_items(
    #[case] model: GtrModelNameCli,
    #[case] items: &[&str],
    #[case] (message, settings): (&str, &[&str]),
  ) {
    let Err(ArgumentError { message: actual_message, settings: actual_settings }) =
      parse_model_params(model, &helpers::owned(items))
    else {
      panic!("expected an argument error");
    };
    assert_eq!((message, helpers::owned(settings)), (actual_message.as_str(), actual_settings));
  }

  #[test]
  fn test_model_params_frequencies_build_an_amino_acid_model() -> Result<(), Report> {
    let pi = Array1::range(1.0, 21.0, 1.0);
    let item = format!("pi={}", pi.iter().join(","));

    let gtr = parse_model_params(GtrModelNameCli::F81, &[item])?.build(AlphabetName::AaNoStop)?;

    pretty_assert_abs_diff_eq!(&pi / pi.sum(), gtr.pi().clone(), epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nuc_too_short(  GtrModelNameCli::HKY85, "pi=0.5,0.5",              AlphabetName::Nuc,      "--model-params item pi: the alphabet nuc has 4 states, so pi needs 4 comma-separated frequencies, one per state, but got 2")]
  #[case::nuc_too_long(   GtrModelNameCli::F81,   "pis=0.2,0.2,0.2,0.2,0.2", AlphabetName::Nuc,      "--model-params item pi: the alphabet nuc has 4 states, so pi needs 4 comma-separated frequencies, one per state, but got 5")]
  #[case::aa_nucleotides( GtrModelNameCli::F81,   "pi=0.1,0.2,0.3,0.4",      AlphabetName::AaNoStop, "--model-params item pi: the alphabet aa-no-stop has 20 states, so pi needs 20 comma-separated frequencies, one per state, but got 4")]
  #[case::tn93_aa(        GtrModelNameCli::TN93,  "Pi=0.1,0.2,0.3,0.4",      AlphabetName::Aa,       "--model-params item pi: the alphabet aa has 21 states, so pi needs 21 comma-separated frequencies, one per state, but got 4")]
  #[trace]
  fn test_model_params_frequencies_need_one_per_state_of_the_alphabet(
    #[case] model: GtrModelNameCli,
    #[case] item: &str,
    #[case] alphabet: AlphabetName,
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let model = parse_model_params(model, &helpers::owned(&[item]))?;
    assert_error!(model.build(alphabet), expected);
    Ok(())
  }

  #[test]
  fn test_model_params_name_models_as_the_config_does() -> Result<(), Report> {
    let expected = GtrModelNameCli::iter()
      .map(|model| format!("\"{model}\""))
      .collect::<Vec<_>>();
    let actual = GtrModelNameCli::iter()
      .map(|model| json_write_str(&model, JsonPretty(false)))
      .collect::<Result<Vec<_>, Report>>()?;
    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::__tests__::test_support::tests::project_root;
    use deser::Deserialize;
    use eyre::Report;
    use ndarray::{Array1, Array2};
    use std::collections::BTreeMap;
    use treetime_utils::adapters::{Array2Rows, ArrayVec};
    use treetime_utils::io::json::json_read_file;

    pub(super) fn owned(items: &[&str]) -> Vec<String> {
      items.iter().map(|item| (*item).to_owned()).collect()
    }

    pub(super) fn v0_models() -> Result<BTreeMap<String, BTreeMap<String, V0Model>>, Report> {
      json_read_file(project_root().join("packages/treetime/src/gtr/__tests__/__fixtures__/gm_gtr_outputs.json"))
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct V0Model {
      pub(super) mu: f64,
      #[deser(deserialize_as = ArrayVec)]
      pub(super) pi: Array1<f64>,
      #[deser(rename = "W", deserialize_as = Array2Rows)]
      pub(super) w: Array2<f64>,
    }
  }
}
