#[cfg(test)]
mod tests {
  use crate::commands::clock::{ClockArgs, run_clock};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::path::PathBuf;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;

  #[test]
  fn test_clock_result_serializes_model_and_regression_results() {
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/ebola/20");
    let outdir = tempfile::tempdir().unwrap();
    let args: ClockArgs = serde_json::from_value(json!({
      "tree": data.join("tree.nwk"),
      "dates": data.join("metadata.tsv"),
      "outdir": outdir.path(),
    }))
    .unwrap();

    let result = run_clock(&args, &NoopCancel, &NoopProgress).unwrap();
    let Value::Object(fields) = serde_json::to_value(&result).unwrap() else {
      panic!("clock result must serialize to a JSON object");
    };

    let shape = fields
      .iter()
      .map(|(key, value)| (key.as_str(), value.is_object(), value.is_array()))
      .collect::<Vec<_>>();
    assert_eq!(
      vec![("clock_model", true, false), ("regression_results", false, true)],
      shape
    );
  }
}
