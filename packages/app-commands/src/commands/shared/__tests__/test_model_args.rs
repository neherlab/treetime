#[cfg(test)]
mod tests {
  use crate::commands::ancestral::aa_model::AaModelName;
  use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
  use crate::commands::shared::model::{GtrModelNameCli, ModelArgs};
  use clap::Parser;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::json;
  use std::iter;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::from_json_value;

  #[rustfmt::skip]
  #[rstest]
  #[case::lower_case(&["--model=jc69",  "--aa-model=jtt92"], (Some(GtrModelNameCli::JC69),  AaModelName::Jtt92))]
  #[case::upper_case(&["--model=JC69",  "--aa-model=JTT92"], (Some(GtrModelNameCli::JC69),  AaModelName::Jtt92))]
  #[case::mixed_case(&["--model=Hky85", "--aa-model=Infer"], (Some(GtrModelNameCli::HKY85), AaModelName::Infer))]
  #[trace]
  fn test_model_args_command_line_ignores_case(
    #[case] args: &[&str],
    #[case] expected: (Option<GtrModelNameCli>, AaModelName),
  ) {
    let raw = TreetimeAncestralArgsRaw::try_parse_from(iter::once(&"ancestral").chain(args)).unwrap();
    assert_eq!(expected, (raw.model_args.model, raw.aa_model));
  }

  #[test]
  fn test_model_args_config_takes_the_exact_value() {
    assert_error!(
      from_json_value::<ModelArgs>(&json!({ "model": "JC69" })),
      "When converting a JSON value: Unexpected: unknown variant `JC69` of GtrModelNameCli, expected one of `infer`, `jc69`, `k80`, `f81`, `hky85`, `t92`, `tn93`, `jtt92`"
    );
  }
}
