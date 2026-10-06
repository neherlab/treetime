#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{GtrModelName, GtrOutput, JC69Params, jc69};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::{JsonPretty, json_write_str};
  use treetime_utils::io::json::{json_read_str, json_value_read_str};
  use treetime_utils::vec_of_owned;

  #[test]
  fn test_gtr_output_without_discrete_states_omits_fields() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let output = GtrOutput::builder().gtr(&gtr).model_name(GtrModelName::JC69).build();
    let json = json_write_str(&output, JsonPretty(false))?;
    let parsed = json_value_read_str(&json)?;

    assert!(
      parsed.get("attribute").is_none(),
      "attribute should be absent when None"
    );
    assert!(parsed.get("states").is_none(), "states should be absent when None");
    assert_eq!(parsed["n_states"], 4);
    assert_eq!(parsed["model_name"], "jc69");
    assert_eq!(parsed["model_type"], "named");

    Ok(())
  }

  #[test]
  fn test_gtr_output_with_discrete_states_includes_fields() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let output = GtrOutput::builder()
      .gtr(&gtr)
      .model_name(GtrModelName::Infer)
      .attribute("country")
      .states(vec_of_owned!["france", "germany", "usa"])
      .build();
    let json = json_write_str(&output, JsonPretty(false))?;
    let parsed = json_value_read_str(&json)?;

    assert_eq!(parsed["attribute"], "country");
    assert_eq!(parsed["states"], serde_json::json!(["france", "germany", "usa"]));
    assert_eq!(parsed["n_states"], 4);
    assert_eq!(parsed["model_type"], "inferred");

    Ok(())
  }

  #[test]
  fn test_gtr_output_with_discrete_states_roundtrip() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let original = GtrOutput::builder()
      .gtr(&gtr)
      .model_name(GtrModelName::Infer)
      .attribute("region")
      .states(vec_of_owned!["asia", "europe"])
      .build();
    let json = json_write_str(&original, JsonPretty(true))?;
    let restored: GtrOutput = json_read_str(&json)?;

    assert_eq!(restored.attribute, Some("region".to_owned()));
    assert_eq!(restored.states, Some(vec!["asia".to_owned(), "europe".to_owned()]));
    assert_eq!(restored.model_type, original.model_type);
    assert_eq!(restored.model_name, original.model_name);
    assert_eq!(restored.n_states, original.n_states);

    Ok(())
  }
}
