#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{
    confidence_by_name, load_gm_mugration_inputs, load_gm_mugration_outputs, run_gm_mugration_case, states_vec,
    trait_assignments_by_name,
  };
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country(          "zika_20_country")]
  #[case::zika_20_country_weights(  "zika_20_country_weights")]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
  #[case::dengue_20_country(        "dengue_20_country")]
  #[case::tb_20_cluster(            "tb_20_cluster")]
  #[case::rsv_a_20_country(         "rsv_a_20_country")]
  #[case::mpox_clade_ii_20_country( "mpox_clade_ii_20_country")]
  #[case::zika_20_missing_values(   "zika_20_country_missing_values")]
  #[case::zika_20_missing_rows(     "zika_20_country_missing_rows")]
  #[trace]
  fn test_gm_mugration_outputs(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs[case];
    let expected = &outputs[case];
    let (output, names) = run_gm_mugration_case(input)?;

    let expected_states = expected.states.clone();
    let actual_states = states_vec(&output);
    assert_eq!(expected_states, actual_states);

    let expected_trait_assignments = expected.trait_assignments.clone();
    let actual_trait_assignments = trait_assignments_by_name(&output, &names);
    assert_eq!(expected_trait_assignments, actual_trait_assignments);

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_20_country(          "zika_20_country")]
  #[case::zika_20_country_weights(  "zika_20_country_weights")]
  #[case::lassa_l_20_country(       "lassa_L_20_country")]
  #[case::dengue_20_country(        "dengue_20_country")]
  #[case::tb_20_cluster(            "tb_20_cluster")]
  #[case::rsv_a_20_country(         "rsv_a_20_country")]
  #[case::mpox_clade_ii_20_country( "mpox_clade_ii_20_country")]
  #[case::zika_20_missing_values(   "zika_20_country_missing_values")]
  #[case::zika_20_missing_rows(     "zika_20_country_missing_rows")]
  #[trace]
  fn test_gm_mugration_confidence_outputs(#[case] case: &str) -> Result<(), Report> {
    let inputs = load_gm_mugration_inputs();
    let outputs = load_gm_mugration_outputs();
    let input = &inputs[case];
    let expected = &outputs[case];
    let (output, names) = run_gm_mugration_case(input)?;

    assert_eq!(expected.states, states_vec(&output));

    let actual_confidence = confidence_by_name(&output, &names);
    for (node_name, expected_profile) in &expected.confidence {
      let actual_profile = actual_confidence
        .get(node_name)
        .unwrap_or_else(|| panic!("missing confidence for node '{node_name}'"));
      let expected_arr = Array1::from_vec(expected_profile.clone());
      pretty_assert_abs_diff_eq!(expected_arr, actual_profile, epsilon = 1e-7);
    }

    Ok(())
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::mugration::pipeline::{MugrationInput, MugrationOutput, MugrationParams, run};
    use crate::mugration::trait_pairing::{TraitPairing, pair_trait_rows};
    use crate::progress::NoopProgress;
    use deser::Deserialize;
    use eyre::Report;
    use indexmap::IndexMap;
    use itertools::Itertools;
    use ndarray::Array1;
    use regex::Regex;
    use std::collections::{BTreeMap, BTreeSet};
    use std::path::PathBuf;
    use std::sync::LazyLock;
    use treetime_graph::graph::Graph;
    use treetime_graph::key_map::NodeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_io::csv::default_name_candidates;
    use treetime_io::discrete_states_csv::discrete_attrs_read_file;
    use treetime_io::nwk::nwk_read_file;
    use treetime_io::state_weights_csv::state_weights_read_file;
    use treetime_utils::io::json::json_read_file;
    use treetime_utils::{make_error, o};

    static GENERATED_NAME: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^NODE_\d{7}$").unwrap());

    #[derive(Debug, Deserialize)]
    pub(super) struct GmMugrationInput {
      tree_path: String,
      metadata_path: String,
      attribute: String,
      name_column: Option<String>,
      parameters: GmMugrationParameters,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct GmMugrationParameters {
      missing_data: String,
      pc: f64,
      sampling_bias_correction: Option<f64>,
      weights_path: Option<String>,
      iterations: usize,
    }

    #[derive(Debug, Deserialize)]
    pub(super) struct GmMugrationOutput {
      pub(crate) states: Vec<String>,
      pub(crate) trait_assignments: BTreeMap<String, String>,
      pub(crate) confidence: BTreeMap<String, Vec<f64>>,
    }

    pub(super) fn load_gm_mugration_inputs() -> IndexMap<String, GmMugrationInput> {
      let path = format!(
        "{}/src/mugration/__tests__/__fixtures__/gm_mugration_inputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      json_read_file(&path).unwrap()
    }

    pub(super) fn load_gm_mugration_outputs() -> IndexMap<String, GmMugrationOutput> {
      let path = format!(
        "{}/src/mugration/__tests__/__fixtures__/gm_mugration_outputs.json",
        env!("CARGO_MANIFEST_DIR")
      );
      json_read_file(&path).unwrap()
    }

    pub(super) fn run_gm_mugration_case(
      fixture: &GmMugrationInput,
    ) -> Result<(MugrationOutput, NodeMap<String>), Report> {
      let project_root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");

      let tree_path = project_root.join(&fixture.tree_path);
      let nwk_parsed = nwk_read_file(&tree_path)?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let metadata_path = project_root.join(&fixture.metadata_path);
      let (rows, _attr_name) = discrete_attrs_read_file::<String>(
        &metadata_path,
        &[',', '\t', ';'],
        &default_name_candidates(),
        fixture.name_column.as_deref(),
        Some(fixture.attribute.as_str()),
        |s| Ok(s.to_owned()),
      )?;
      let TraitPairing {
        traits,
        observed_values,
        ..
      } = pair_trait_rows(rows, &graph, &names);

      let weights = match &fixture.parameters.weights_path {
        Some(weights_path) => {
          let weights_filepath = project_root.join(weights_path);
          Some(
            state_weights_read_file(&weights_filepath, &[',', '\t', ';'])?
              .into_iter()
              .collect(),
          )
        },
        None => None,
      };

      let params = MugrationParams {
        missing_data: fixture.parameters.missing_data.clone(),
        pc: fixture.parameters.pc,
        missing_weights_threshold: 0.5,
        iterations: fixture.parameters.iterations,
        sampling_bias_correction: fixture.parameters.sampling_bias_correction,
        smooth_initial_pi: false,
        filter_uninformative_root: false,
      };
      let input = MugrationInput {
        graph,
        traits,
        traits_source: o!("metadata"),
        observed_values,
        weights,
        branch_lengths,
      };
      let output = run(&params, input, &names, &NoopCancel, &NoopProgress).map_err(|err| err.into_report())?;
      let v0_names = v0_node_names(&output.graph, &names)?;
      Ok((output, v0_names))
    }

    pub(super) fn states_vec(output: &MugrationOutput) -> Vec<String> {
      output.states.iter().map(|s| s.to_owned()).collect()
    }

    pub(super) fn trait_assignments_by_name(
      output: &MugrationOutput,
      names: &NodeMap<String>,
    ) -> BTreeMap<String, String> {
      names
        .iter()
        .map(|(key, name)| {
          let state = output.states.name(output.posteriors[key].state()).unwrap();
          (name.clone(), state.to_owned())
        })
        .collect()
    }

    pub(super) fn confidence_by_name(
      output: &MugrationOutput,
      names: &NodeMap<String>,
    ) -> BTreeMap<String, Array1<f64>> {
      names
        .iter()
        .map(|(key, name)| (name.clone(), output.posteriors[key].profile().clone()))
        .collect()
    }

    fn v0_node_names(graph: &Graph, names: &NodeNames) -> Result<NodeMap<String>, Report> {
      let mut v0_names = NodeMap::new();
      for node in graph.get_nodes().filter(|node| node.is_leaf()) {
        v0_names.insert(node.key(), names[node.key()].to_owned());
      }
      for node in graph.get_nodes().filter(|node| !node.is_leaf()) {
        let name = &names[node.key()];
        if !GENERATED_NAME.is_match(name) {
          return make_error!(
            "Internal node label '{name}': the v0 node names of the harness need a tree without internal labels"
          );
        }
      }

      let mut used: BTreeSet<String> = v0_names.values().cloned().collect();
      let mut order = vec![];
      ladderized_internal_preorder(graph, graph.get_exactly_one_root()?.key(), &mut order);
      let mut internal_node_count = 0;
      for key in order {
        let mut name = format!("NODE_{internal_node_count:07}");
        while used.contains(&name) {
          internal_node_count += 1;
          name = format!("NODE_{internal_node_count:07}");
        }
        used.insert(name.clone());
        v0_names.insert(key, name);
        internal_node_count += 1;
      }
      Ok(v0_names)
    }

    fn ladderized_internal_preorder(graph: &Graph, key: GraphNodeKey, order: &mut Vec<GraphNodeKey>) {
      let node = graph.get_node(key).unwrap();
      if node.is_leaf() {
        return;
      }
      order.push(key);
      let children = graph
        .children_keys_of(node)
        .map(|(child, _)| child)
        .sorted_by_key(|&child| leaf_count(graph, child));
      for child in children {
        ladderized_internal_preorder(graph, child, order);
      }
    }

    fn leaf_count(graph: &Graph, key: GraphNodeKey) -> usize {
      let node = graph.get_node(key).unwrap();
      if node.is_leaf() {
        1
      } else {
        graph
          .children_keys_of(node)
          .map(|(child, _)| leaf_count(graph, child))
          .sum()
      }
    }
  }
}
