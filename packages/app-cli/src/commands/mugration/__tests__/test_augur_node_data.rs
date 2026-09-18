#[cfg(test)]
mod tests {
  use app_output::augur_node_data_mugration::build_augur_node_data_json;
  use app_output::mugration_result::MugrationResult;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime::cancel::NoopCancel;
  use treetime::mugration::pipeline::{self, MugrationInput, MugrationParams};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
  use treetime_utils::o;
  use util_augur_node_data_json::AugurNodeDataJsonTraits;

  fn run_and_serialize(tree: &str, traits: &std::collections::BTreeMap<String, String>) -> String {
    let nwk_parsed = nwk_read_str(tree).unwrap();
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let params = MugrationParams {
      missing_data: o!("?"),
      pc: None,
      missing_weights_threshold: 0.5,
      iterations: 5,
      sampling_bias_correction: None,
      smooth_initial_pi: false,
      filter_uninformative_root: false,
    };
    let input = MugrationInput {
      graph,
      traits: traits.clone(),
      weights: None,
      branch_lengths: branch_lengths.clone(),
    };
    let output = pipeline::run(&params, input, &names, &NoopCancel).unwrap();
    let result = MugrationResult::new(&output, &confidences, &names, &branch_lengths, "country");
    let data = build_augur_node_data_json(&result, &output).unwrap();
    json_write_str(&data, JsonPretty(true)).unwrap()
  }

  #[test]
  fn test_augur_node_data_mugration_full_output() {
    let actual = run_and_serialize(
      "(A:0.1,B:0.2)root;",
      &btreemap! { o!("A") => o!("usa"), o!("B") => o!("germany") },
    );

    let expected = format!(
      r#"{{
  "generated_by": {{
    "program": "treetime",
    "version": "{version}"
  }},
  "nodes": {{
    "A": {{
      "country": "usa",
      "country_confidence": {{
        "usa": 1.0
      }},
      "country_entropy": -1.000088900581841e-12
    }},
    "B": {{
      "country": "germany",
      "country_confidence": {{
        "germany": 1.0
      }},
      "country_entropy": -1.000088900581841e-12
    }},
    "root": {{
      "country": "usa",
      "country_confidence": {{
        "germany": 0.31787306601445664,
        "usa": 0.6821269339855435
      }},
      "country_entropy": 0.6252558275846285
    }}
  }},
  "models": {{
    "country": {{
      "rate": 1.9605478705568704,
      "alphabet": [
        "germany",
        "usa",
        "?"
      ],
      "equilibrium_probabilities": [
        0.40794069765038415,
        0.5920593023496159
      ],
      "transition_matrix": [
        [
          0.0,
          2.0701783432087786
        ],
        [
          2.0701783432087786,
          0.0
        ]
      ]
    }}
  }},
  "branches": {{
    "B": {{
      "labels": {{
        "country": "usa → germany"
      }}
    }},
    "root": {{
      "labels": {{
        "country": "usa"
      }}
    }}
  }}
}}"#,
      version = env!("CARGO_PKG_VERSION")
    );

    assert_eq!(expected, actual.trim());
  }

  #[test]
  fn test_augur_node_data_mugration_roundtrip() {
    let json_str = run_and_serialize(
      "(A:0.1,B:0.2)root;",
      &btreemap! { o!("A") => o!("usa"), o!("B") => o!("germany") },
    );
    let original: serde_json::Value = serde_json::from_str(&json_str).unwrap();
    let typed: AugurNodeDataJsonTraits = json_read_str(&json_str).unwrap();
    let roundtripped: serde_json::Value = serde_json::to_value(&typed).unwrap();
    assert_eq!(original, roundtripped);
  }
}
