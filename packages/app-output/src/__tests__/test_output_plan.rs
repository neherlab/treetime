#[cfg(test)]
mod tests {
  use crate::output_plan::OutputSelection;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::str::FromStr;
  use strum::IntoEnumIterator;

  #[test]
  fn test_output_plan_selection_tag_matches_serde_name() {
    let serde_names = OutputSelection::iter()
      .map(|selection| serde_json::to_value(selection).unwrap())
      .collect_vec();
    let tag_names = OutputSelection::iter()
      .map(|selection| serde_json::Value::String(selection.as_ref().to_owned()))
      .collect_vec();
    assert_eq!(serde_names, tag_names);
  }

  #[test]
  fn test_output_plan_selection_from_str_roundtrip() {
    let expected = OutputSelection::iter().collect_vec();
    let actual = OutputSelection::iter()
      .map(|selection| OutputSelection::from_str(selection.as_ref()).unwrap())
      .collect_vec();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_output_plan_selection_from_str_rejects_unknown_name() {
    assert_eq!(
      Err(strum::ParseError::VariantNotFound),
      OutputSelection::from_str("mat_pb")
    );
  }
}
