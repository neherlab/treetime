#[cfg(test)]
mod tests {
  use crate::ancestral::params::AncestralParams;
  use crate::ancestral::params::MethodAncestral;
  use crate::ancestral::plan::{ReconPlan, resolve_plan};
  use crate::error::OperationError;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::create::Representation;
  use crate::partition::marginal::sample::SampleMode;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[test]
  fn test_plan_parsimony_resolves_to_fitch() {
    let params = helpers::params(MethodAncestral::Parsimony);
    assert!(matches!(resolve_plan(&params), Ok(ReconPlan::Fitch)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default_is_sparse(None,        false)]
  #[case::dense(            Some(true),  true)]
  #[case::sparse(           Some(false), false)]
  #[trace]
  fn test_plan_marginal_representation(#[case] dense: Option<bool>, #[case] expected_dense: bool) {
    let params = AncestralParams {
      dense,
      ..helpers::params(MethodAncestral::Marginal)
    };
    let Ok(ReconPlan::Marginal { representation, .. }) = resolve_plan(&params) else {
      panic!("marginal method must resolve to a marginal plan");
    };
    assert_eq!(expected_dense, matches!(representation, Representation::Dense));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::infer_with_iterations(3, GtrModelName::Infer, Some(3))]
  #[case::infer_without_iterations(0, GtrModelName::Infer, None)]
  #[case::named_model_with_iterations(3, GtrModelName::JC69, None)]
  #[trace]
  fn test_plan_marginal_gtr_refinement_gate(
    #[case] gtr_iterations: usize,
    #[case] model: GtrModelName,
    #[case] expected: Option<usize>,
  ) {
    let params = AncestralParams {
      model,
      gtr_iterations,
      ..helpers::params(MethodAncestral::Marginal)
    };
    let Ok(ReconPlan::Marginal { gtr_refinement, .. }) = resolve_plan(&params) else {
      panic!("marginal method must resolve to a marginal plan");
    };
    assert_eq!(expected, gtr_refinement);
  }

  #[test]
  fn test_plan_marginal_keeps_model() {
    let params = AncestralParams {
      model: GtrModelName::HKY85,
      ..helpers::params(MethodAncestral::Marginal)
    };
    let Ok(ReconPlan::Marginal { model, .. }) = resolve_plan(&params) else {
      panic!("marginal method must resolve to a marginal plan");
    };
    assert!(matches!(model, GtrModelName::HKY85));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::site_specific_gtr(
    MethodAncestral::Marginal, true, SampleMode::Argmax,
    "--site-specific-gtr is not implemented"
  )]
  #[case::site_specific_gtr_before_joint(
    MethodAncestral::Joint, true, SampleMode::Root,
    "--site-specific-gtr is not implemented"
  )]
  #[case::sampling_with_parsimony(
    MethodAncestral::Parsimony, false, SampleMode::Root,
    "--sample-from-profile=Root requires --method-anc=marginal. Posterior sampling is only defined for marginal reconstruction; Parsimony has no posterior profile to sample."
  )]
  #[case::sampling_before_joint(
    MethodAncestral::Joint, false, SampleMode::All,
    "--sample-from-profile=All requires --method-anc=marginal. Posterior sampling is only defined for marginal reconstruction; Joint has no posterior profile to sample."
  )]
  #[case::joint(
    MethodAncestral::Joint, false, SampleMode::Argmax,
    "Joint ancestral reconstruction has been removed. Available methods: marginal, parsimony"
  )]
  #[trace]
  fn test_plan_rejects_invalid_params(
    #[case] method: MethodAncestral,
    #[case] site_specific_gtr: bool,
    #[case] sample_from_profile: SampleMode,
    #[case] expected: &str,
  ) {
    let params = AncestralParams {
      site_specific_gtr,
      sample_from_profile,
      ..helpers::params(method)
    };
    let Err(OperationError::InvalidParams(report)) = resolve_plan(&params) else {
      panic!("expected an invalid-params rejection");
    };
    assert_eq!(expected, report.to_string());
  }

  mod helpers {
    use super::*;

    pub(super) fn params(method: MethodAncestral) -> AncestralParams {
      AncestralParams {
        method,
        model: GtrModelName::Infer,
        dense: None,
        include_leaves: false,
        impute_missing_data: false,
        gtr_iterations: 0,
        site_specific_gtr: false,
        seed: None,
        sample_from_profile: SampleMode::Argmax,
      }
    }
  }
}
