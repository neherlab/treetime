#[cfg(test)]
mod tests {
  use crate::commands::mugration::discrete_marginal::attach_traits;
  use crate::o;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_io::json::{JsonPretty, json_write_str};
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_discrete_marginal_attach_traits_maps_observed_and_missing_profiles() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("?"),
    };

    attach_traits(&mut partition, &graph, &traits)?;

    let node_a = helpers::get_node_data(&graph, &partition, "A");
    let actual_node_a = json_write_str(node_a, JsonPretty(true))?;
    let expected_node_a = o!(indoc! {r#"{
      "observed": 1,
      "profile": [
        0.0,
        1.0
      ],
      "log_lh": 0.0
    }"#});
    assert_eq!(expected_node_a, actual_node_a);

    let node_b = helpers::get_node_data(&graph, &partition, "B");
    let actual_node_b = json_write_str(node_b, JsonPretty(true))?;
    let expected_node_b = o!(indoc! {r#"{
      "observed": null,
      "profile": [
        0.5,
        0.5
      ],
      "log_lh": 0.0
    }"#});
    assert_eq!(expected_node_b, actual_node_b);

    assert_eq!(graph.get_edges().len(), partition.edges.len());

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_rejects_tree_leaf_missing_from_metadata() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
    };

    let error = attach_traits(&mut partition, &graph, &traits).unwrap_err().to_string();

    assert!(error.contains("tree leaves missing from metadata"), "{error}");
    assert!(error.contains('B'), "{error}");

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_rejects_metadata_name_missing_from_tree() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = btreemap! {
      o!("A") => o!("usa"),
      o!("B") => o!("germany"),
      o!("C") => o!("usa"),
    };

    let error = attach_traits(&mut partition, &graph, &traits).unwrap_err().to_string();

    assert!(error.contains("metadata names missing from tree leaves"), "{error}");
    assert!(error.contains('C'), "{error}");

    Ok(())
  }

  mod helpers {
    use crate::gtr::gtr::{GTR, GTRParams};
    use crate::representation::discrete_states::DiscreteStates;
    use crate::representation::partition::discrete::PartitionDiscrete;
    use crate::representation::payload::ancestral::GraphAncestral;
    use crate::representation::payload::discrete::DiscreteNodeData;
    use crate::test_utils::find_node_key_by_name;
    use eyre::Report;
    use ndarray::Array1;

    pub(super) fn make_partition(states: [&str; 2]) -> Result<PartitionDiscrete, Report> {
      let discrete_states = DiscreteStates::from_values(states.into_iter(), "?");
      let n_states = discrete_states.len();
      let gtr = GTR::new(GTRParams {
        n_states,
        mu: 1.0,
        W: None,
        pi: Array1::from_elem(n_states, 1.0 / n_states as f64),
      })?;

      Ok(PartitionDiscrete::new(0, gtr, discrete_states))
    }

    pub(super) fn get_node_data<'a>(
      graph: &GraphAncestral,
      partition: &'a PartitionDiscrete,
      name: &str,
    ) -> &'a DiscreteNodeData {
      let node_key = find_node_key_by_name(graph, name).unwrap_or_else(|| panic!("Missing test node '{name}'"));
      partition
        .nodes
        .get(&node_key)
        .unwrap_or_else(|| panic!("Missing discrete node data for '{name}'"))
    }
  }
}
