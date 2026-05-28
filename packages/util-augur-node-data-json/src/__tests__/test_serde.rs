#[cfg(test)]
mod tests {
  use crate::*;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  #[test]
  fn test_ancestral_serialization_minimal() {
    let data = AugurNodeDataJsonAncestral {
      generated_by: Some(AugurNodeDataJsonGeneratedBy {
        program: "treetime".to_owned(),
        version: "1.0.0".to_owned(),
      }),
      metadata: AugurNodeDataJsonAncestralMeta {
        annotations: Some(AugurNodeDataJsonAnnotations {
          nuc: Some(AugurNodeDataJsonAnnotationEntry {
            start: Some(1),
            end: Some(50),
            strand: Some("+".to_owned()),
            entry_type: Some("source".to_owned()),
            segments: None,
            other: btreemap! {},
          }),
          other: btreemap! {},
        }),
        reference: Some(btreemap! {
          "nuc".to_owned() => "ACGT".to_owned(),
        }),
        mask: Some("0000".to_owned()),
        other: btreemap! {},
      },
      nodes: btreemap! {
        "root".to_owned() => AugurNodeDataJsonAncestralNode {
          muts: vec![],
          sequence: Some("ACGT".to_owned()),
          aa_muts: None,
          aa_sequences: None,
          other: btreemap! {},
        },
        "tip_A".to_owned() => AugurNodeDataJsonAncestralNode {
          muts: vec!["A1T".to_owned(), "C3G".to_owned()],
          sequence: Some("TCGG".to_owned()),
          aa_muts: None,
          aa_sequences: None,
          other: btreemap! {},
        },
      },
    };

    let json = json_write_str(&data, JsonPretty(true)).unwrap();
    assert!(json.contains("\"generated_by\""));
    assert!(json.contains("\"nodes\""));
    assert!(json.contains("\"annotations\""));
    assert!(json.contains("\"mask\": \"0000\""));
    assert!(json.contains("\"reference\""));
    assert!(json.contains("\"muts\": []"));
    assert!(json.contains("\"sequence\": \"ACGT\""));

    let roundtrip: AugurNodeDataJsonAncestral = serde_json::from_str(&json).unwrap();
    assert_eq!(2, roundtrip.nodes.len());
    assert_eq!(Vec::<String>::new(), roundtrip.nodes["root"].muts);
    assert_eq!(Some("ACGT".to_owned()), roundtrip.nodes["root"].sequence);
    assert_eq!(vec!["A1T", "C3G"], roundtrip.nodes["tip_A"].muts);
    assert_eq!(Some("0000".to_owned()), roundtrip.metadata.mask);
  }

  #[test]
  fn test_refine_serialization_clock() {
    let data = AugurNodeDataJsonRefine {
      generated_by: Some(AugurNodeDataJsonGeneratedBy {
        program: "treetime".to_owned(),
        version: "1.0.0".to_owned(),
      }),
      metadata: AugurNodeDataJsonRefineMeta {
        alignment: None,
        input_tree: None,
        clock: Some(AugurNodeDataJsonClock {
          rate: 0.003,
          intercept: -2.21,
          rtt_tmrca: 736.6666666666666,
          cov: Some(vec![vec![1e-8, 0.0], vec![0.0, 0.5]]),
          rate_std: Some(1e-4),
          other: btreemap! {},
        }),
        other: btreemap! {},
      },
      nodes: btreemap! {
        "leaf1".to_owned() => AugurNodeDataJsonRefineNode {
          branch_length: 0.001,
          numdate: Some(2013.999),
          date: Some("2013-12-31".to_owned()),
          raw_date: Some("2013-12-XX".to_owned()),
          date_inferred: Some(true),
          num_date_confidence: Some([2013.5, 2014.5]),
          ..Default::default()
        },
      },
    };

    let json = json_write_str(&data, JsonPretty(true)).unwrap();
    assert!(json.contains("\"rtt_Tmrca\""));
    assert!(json.contains("\"numdate\""));
    assert!(json.contains("\"raw_date\""));
    assert!(json.contains("\"num_date_confidence\""));
  }

  #[test]
  fn test_traits_serialization() {
    let data = AugurNodeDataJsonTraits {
      generated_by: Some(AugurNodeDataJsonGeneratedBy {
        program: "treetime".to_owned(),
        version: "1.0.0".to_owned(),
      }),
      metadata: AugurNodeDataJsonTraitsMeta {
        models: Some(btreemap! {
          "region".to_owned() => AugurNodeDataJsonTraitModel {
            rate: 3.737,
            alphabet: vec!["North America".to_owned(), "Oceania".to_owned(), "?".to_owned()],
            equilibrium_probabilities: vec![0.5, 0.5],
            transition_matrix: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            other: btreemap! {},
          },
        }),
        other: btreemap! {},
      },
      nodes: btreemap! {
        "root".to_owned() => AugurNodeDataJsonTraitsNode {
          fields: btreemap! {
            "region".to_owned() => serde_json::Value::String("North America".to_owned()),
            "region_confidence".to_owned() => serde_json::json!({"North America": 0.95, "Oceania": 0.05}),
            "region_entropy".to_owned() => serde_json::json!(0.286),
          },
        },
      },
    };

    let json = json_write_str(&data, JsonPretty(true)).unwrap();
    assert!(json.contains("\"models\""));
    assert!(json.contains("\"region_confidence\""));
    assert!(json.contains("\"equilibrium_probabilities\""));
  }

  #[test]
  fn test_roundtrip_ancestral() {
    let data = AugurNodeDataJsonAncestral {
      generated_by: Some(AugurNodeDataJsonGeneratedBy {
        program: "augur".to_owned(),
        version: "25.0.0".to_owned(),
      }),
      metadata: AugurNodeDataJsonAncestralMeta {
        annotations: None,
        reference: None,
        mask: None,
        other: btreemap! {},
      },
      nodes: btreemap! {
        "node1".to_owned() => AugurNodeDataJsonAncestralNode {
          muts: vec!["A5C".to_owned()],
          sequence: None,
          aa_muts: None,
          aa_sequences: None,
          other: btreemap! {},
        },
      },
    };

    let json = json_write_str(&data, JsonPretty(false)).unwrap();
    let roundtrip: AugurNodeDataJsonAncestral = serde_json::from_str(&json).unwrap();

    assert_eq!(data.generated_by.unwrap().program, roundtrip.generated_by.unwrap().program);
    assert_eq!(data.nodes["node1"].muts, roundtrip.nodes["node1"].muts);
  }

  #[test]
  fn test_forward_compat_unknown_fields() {
    let json = r#"{
      "generated_by": {"program": "augur", "version": "99.0.0"},
      "nodes": {"n1": {"muts": [], "extra_field": 42}},
      "annotations": {"nuc": {"start": 1, "end": 100, "strand": "+", "type": "source"}},
      "future_top_level_key": "preserved"
    }"#;

    let parsed: AugurNodeDataJsonAncestral = serde_json::from_str(json).unwrap();
    assert_eq!(parsed.nodes["n1"].other["extra_field"], serde_json::json!(42));
    assert_eq!(parsed.metadata.other["future_top_level_key"], serde_json::json!("preserved"));
  }

  #[test]
  fn test_skip_none_fields() {
    let node = AugurNodeDataJsonAncestralNode {
      muts: vec![],
      sequence: None,
      aa_muts: None,
      aa_sequences: None,
      other: btreemap! {},
    };

    let json = json_write_str(&node, JsonPretty(false)).unwrap();
    assert!(!json.contains("sequence"));
    assert!(!json.contains("aa_muts"));
    assert!(!json.contains("aa_sequences"));
  }
}
