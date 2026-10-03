#[cfg(test)]
mod tests {
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use treetime::clock::clock_model::ClockModel;
  use treetime::gtr::get_gtr::{GtrModelName, GtrOutput, JC69Params, jc69};
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};

  #[test]
  fn test_model_json_gtr_output_round_trips() -> Result<(), Report> {
    let gtr = jc69(JC69Params::default())?;
    let output = GtrOutput::builder().gtr(&gtr).model_name(GtrModelName::JC69).build();

    let written = json_write_str(&output, JsonPretty(true))?;
    let read: GtrOutput = json_read_str(&written)?;

    assert_eq!(written, json_write_str(&read, JsonPretty(true))?);
    let fields: Vec<String> = json_read_str::<Value>(&written)?
      .as_object()
      .expect("GTR JSON is an object")
      .keys()
      .cloned()
      .collect();
    assert_eq!(vec!["model_type", "model_name", "mu", "pi", "W", "n_states"], fields);
    Ok(())
  }

  #[test]
  fn test_model_json_clock_model_round_trips() -> Result<(), Report> {
    let expected = json!({
      "clock_rate": 0.001,
      "intercept": -2.0,
      "stats": {
        "estimated": {
          "chisq": 0.5,
          "r_val": 0.9,
          "hessian": [[4.0, 2.0], [2.0, 3.0]],
          "cov": [[0.375, -0.25], [-0.25, 0.5]],
        },
      },
    });

    let read: ClockModel = json_read_str(expected.to_string())?;
    let written: Value = json_read_str(&json_write_str(&read, JsonPretty(true))?)?;

    assert_eq!(expected, written);
    Ok(())
  }
}
