#[cfg(test)]
mod tests {
  use crate::config::resolve_paths::{resolve_config_paths, resolve_config_paths_where};
  use helpers::schema;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use std::path::Path;

  const BASE: &str = "/data/zika";

  #[rustfmt::skip]
  #[rstest]
  #[case::input(          json!({ "tree": "tree.nwk" }),                   json!({ "tree": "/data/zika/tree.nwk" }))]
  #[case::input_template( json!({ "translations": "aa/{cds}.fasta" }),     json!({ "translations": "/data/zika/aa/{cds}.fasta" }))]
  #[case::output(         json!({ "output_all": "out" }),                  json!({ "output_all": "/data/zika/out" }))]
  #[case::parent_folder(  json!({ "tree": "../shared/tree.nwk" }),         json!({ "tree": "/data/zika/../shared/tree.nwk" }))]
  #[case::absolute(       json!({ "tree": "/elsewhere/tree.nwk" }),        json!({ "tree": "/elsewhere/tree.nwk" }))]
  #[case::stdin_in_array( json!({ "alignment": ["-", "aln.fasta"] }),      json!({ "alignment": ["-", "/data/zika/aln.fasta"] }))]
  #[case::stdout_output(  json!({ "output_all": "-" }),                    json!({ "output_all": "-" }))]
  #[case::empty_string(   json!({ "tree": "" }),                           json!({ "tree": "" }))]
  #[case::tilde(          json!({ "tree": "~/tree.nwk" }),                 json!({ "tree": "/data/zika/~/tree.nwk" }))]
  #[case::not_a_string(   json!({ "tree": 7 }),                            json!({ "tree": 7 }))]
  #[case::not_a_path(     json!({ "attribute": "country.tsv" }),           json!({ "attribute": "country.tsv" }))]
  #[trace]
  fn test_resolve_paths_rewrites_relative_paths_of_path_settings(#[case] config: Value, #[case] expected: Value) {
    let mut config = config;
    resolve_config_paths(&mut config, &schema(), Path::new(BASE)).unwrap();
    assert_eq!(expected, config);
  }

  #[test]
  fn test_resolve_paths_leaves_a_document_that_is_not_a_mapping() {
    let mut config = json!(["tree.nwk"]);
    resolve_config_paths(&mut config, &schema(), Path::new(BASE)).unwrap();
    assert_eq!(json!(["tree.nwk"]), config);
  }

  #[test]
  fn test_resolve_paths_where_skips_settings_outside_the_selection() {
    let mut config = json!({ "tree": "tree.nwk", "alignment": ["aln.fasta"] });
    resolve_config_paths_where(&mut config, &schema(), Path::new(BASE), |key| key == "tree").unwrap();
    assert_eq!(
      json!({ "tree": "/data/zika/tree.nwk", "alignment": ["aln.fasta"] }),
      config
    );
  }

  mod helpers {
    use serde_json::{Value, json};

    pub(super) fn schema() -> Value {
      json!({
        "type": "object",
        "properties": {
          "tree": { "type": "string", "x-path": "input" },
          "alignment": { "type": "array", "items": { "type": "string" }, "x-path": "input" },
          "translations": { "type": "string", "x-path": "input-template" },
          "output_all": { "type": "string", "x-path": "output" },
          "attribute": { "type": "string" }
        }
      })
    }
  }
}
