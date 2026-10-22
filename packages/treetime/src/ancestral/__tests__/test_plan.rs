#[cfg(test)]
mod tests {
  use crate::ancestral::params::AncestralParams;
  use crate::ancestral::params::MethodAncestral;
  use crate::ancestral::plan::{ReconstructionPlan, resolve_plan};
  use crate::error::OperationError;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::create::Representation;
  use crate::partition::marginal::sample::SampleMode;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[test]
  fn test_plan_parsimony_resolves_to_fitch() {
    let params = helpers::params(MethodAncestral::Parsimony);
    assert!(matches!(resolve_plan(&params), Ok(ReconstructionPlan::Fitch)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unset_is_sparse_known_issue_n_representation_infer_dense_stub(None, false)]
  #[case::dense(Some(true), true)]
  #[case::sparse(Some(false), false)]
  #[trace]
  fn test_plan_marginal_representation(#[case] dense: Option<bool>, #[case] expected_dense: bool) {
    let params = AncestralParams {
      dense,
      ..helpers::params(MethodAncestral::Marginal)
    };
    let Ok(ReconstructionPlan::Marginal { representation, .. }) = resolve_plan(&params) else {
      panic!("marginal method must resolve to a marginal plan");
    };
    assert_eq!(expected_dense, matches!(representation, Representation::Dense));
  }

  #[test]
  fn test_plan_marginal_keeps_model() {
    let params = AncestralParams {
      model: GtrModelName::HKY85,
      ..helpers::params(MethodAncestral::Marginal)
    };
    let Ok(ReconstructionPlan::Marginal { model, .. }) = resolve_plan(&params) else {
      panic!("marginal method must resolve to a marginal plan");
    };
    assert!(matches!(model, GtrModelName::HKY85));
  }

  #[test]
  fn test_plan_rejects_sampling_with_parsimony() {
    let params = AncestralParams {
      sample_from_profile: SampleMode::Root,
      ..helpers::params(MethodAncestral::Parsimony)
    };
    let Err(OperationError::InvalidParams(report)) = resolve_plan(&params) else {
      panic!("expected an invalid-params rejection");
    };
    assert_eq!(
      "--sample-from-profile=Root requires --method-anc=marginal. Posterior sampling is only defined for marginal reconstruction; Parsimony has no posterior profile to sample.",
      report.to_string()
    );
  }

  mod helpers {
    use super::*;

    pub(super) fn params(method: MethodAncestral) -> AncestralParams {
      AncestralParams {
        method,
        model: GtrModelName::Infer,
        dense: None,
        include_leaves: false,
        report_ambiguous: true,
        impute_missing_data: false,
        seed: 0,
        sample_from_profile: SampleMode::Argmax,
      }
    }
  }
}
