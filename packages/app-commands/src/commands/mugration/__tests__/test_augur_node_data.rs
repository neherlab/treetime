#[cfg(test)]
mod tests {
  use app_output::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, Divergence, TreeTraits};
  use app_output::augur_node_data_traits::build_augur_node_data_traits;
  use helpers::run_and_serialize;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use serde_json::Value;
  use treetime::cancel::NoopCancel;
  use treetime::gtr::infer_gtr::common::DEFAULT_PSEUDO_COUNT;
  use treetime::mugration::pipeline::{self, MugrationInput, MugrationParams};
  use treetime::mugration::trait_pairing::{TraitPairing, pair_trait_rows};
  use treetime::progress::NoopProgress;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::io::json::{JsonPretty, json_read_str, json_value_read_str, json_write_str, to_json_value};
  use treetime_utils::{o, pretty_assert_ulps_eq};
  use util_augur_node_data_json::AugurNodeDataJsonTraits;

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
        "germany": 0.3219171495641157,
        "usa": 0.6780828504358842
      }},
      "country_entropy": 0.6283061142879484
    }}
  }},
  "models": {{
    "country": {{
      "rate": 1.5295225032725375,
      "alphabet": [
        "germany",
        "usa",
        "?"
      ],
      "equilibrium_probabilities": [
        0.41863862983617267,
        0.5813613701638274
      ],
      "transition_matrix": [
        [
          0.0,
          2.054397761926173
        ],
        [
          2.054397761926173,
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
  fn test_augur_node_data_mugration_model_and_root_match_v0() {
    let actual = json_value_read_str(run_and_serialize(
      "(A:0.1,B:0.2)root;",
      &btreemap! { o!("A") => o!("usa"), o!("B") => o!("germany") },
    ))
    .unwrap();
    let model = &actual["models"]["country"];
    let root = &actual["nodes"]["root"]["country_confidence"];
    let number = |value: &Value| value.as_f64().unwrap();

    pretty_assert_ulps_eq!(1.5295225032725375, number(&model["rate"]), max_ulps = 2);
    pretty_assert_ulps_eq!(
      0.4186386298361726,
      number(&model["equilibrium_probabilities"][0]),
      max_ulps = 2
    );
    pretty_assert_ulps_eq!(
      0.5813613701638275,
      number(&model["equilibrium_probabilities"][1]),
      max_ulps = 2
    );
    pretty_assert_ulps_eq!(
      2.0543977619261726,
      number(&model["transition_matrix"][0][1]),
      max_ulps = 2
    );
    pretty_assert_ulps_eq!(0.3219171495641157, number(&root["germany"]), max_ulps = 2);
    pretty_assert_ulps_eq!(0.6780828504358842, number(&root["usa"]), max_ulps = 2);
  }

  #[test]
  fn test_augur_node_data_mugration_roundtrip() {
    let json_str = run_and_serialize(
      "(A:0.1,B:0.2)root;",
      &btreemap! { o!("A") => o!("usa"), o!("B") => o!("germany") },
    );
    let original = json_value_read_str(&json_str).unwrap();
    let typed: AugurNodeDataJsonTraits = json_read_str(&json_str).unwrap();
    let roundtripped = to_json_value(&typed).unwrap();
    assert_eq!(original, roundtripped);
  }

  mod helpers {
    use super::*;

    pub(super) fn run_and_serialize(tree: &str, traits: &std::collections::BTreeMap<String, String>) -> String {
      let nwk_parsed = nwk_read(tree.as_bytes()).unwrap();
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let params = MugrationParams {
        missing_data: o!("?"),
        pc: DEFAULT_PSEUDO_COUNT,
        missing_weights_threshold: 0.5,
        iterations: 5,
        sampling_bias_correction: None,
        smooth_initial_pi: false,
        filter_uninformative_root: false,
      };
      let TraitPairing {
        traits,
        observed_values,
        ..
      } = pair_trait_rows(traits.clone(), &graph, &names);
      let input = MugrationInput {
        graph,
        traits,
        traits_source: o!("metadata 'metadata.tsv'"),
        observed_values,
        weights: None,
        branch_lengths,
      };
      let output = pipeline::run(&params, input, &names, &NoopCancel, &NoopProgress).unwrap();
      let branch_lengths = output.graph.get_edges().map(|edge| (edge.key(), None)).collect();
      let annotated = AnnotatedGraph {
        graph: &output.graph,
        names: &names,
        divergence_branch_lengths: &branch_lengths,
        time_branch_lengths: None,
        divergence: Divergence::CumulativeBranchLength,
        sequences: None,
        dates: None,
        traits: Some(TreeTraits {
          attribute: "country",
          states: &output.states,
          posteriors: &output.posteriors,
        }),
      };
      let data = build_augur_node_data_traits(&AnnotatedTreeView::new(&annotated).unwrap(), &output.gtr).unwrap();
      json_write_str(&data, JsonPretty(true)).unwrap()
    }
  }
}
