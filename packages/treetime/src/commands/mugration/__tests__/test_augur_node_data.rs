#[cfg(test)]
mod tests {
  use crate::commands::mugration::augur_node_data::build_augur_node_data_json;
  use crate::mugration::mugration::execute_mugration;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
  use treetime_utils::o;
  use util_augur_node_data_json::AugurNodeDataJsonTraits;

  fn run_and_serialize(tree: &str, traits: &std::collections::BTreeMap<String, String>) -> String {
    let graph = nwk_read_str(tree).unwrap();
    let result = execute_mugration(graph, traits, "country", None, "?", None, 0.5, 5, None, false, false).unwrap();
    let data = build_augur_node_data_json(&result).unwrap();
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
        "germany": 0.31787306595142145,
        "usa": 0.6821269340485786
      }},
      "country_entropy": 0.6252558275364971
    }}
  }},
  "models": {{
    "country": {{
      "rate": 1.9605478625169048,
      "alphabet": [
        "germany",
        "usa",
        "?"
      ],
      "equilibrium_probabilities": [
        0.40794069766081037,
        0.5920593023391897
      ],
      "transition_matrix": [
        [
          0.0,
          2.0701783431923246
        ],
        [
          2.0701783431923246,
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
