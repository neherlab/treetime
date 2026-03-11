#[cfg(test)]
mod tests {
  use crate::commands::mugration::discrete_marginal::attach_traits;
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::Array1;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_discrete_marginal_attach_traits_maps_observed_and_missing_profiles() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = BTreeMap::from([("A".to_owned(), "usa".to_owned()), ("B".to_owned(), "?".to_owned())]);

    attach_traits(&mut partition, &graph, &traits)?;

    let node_a = helpers::get_node_data(&graph, &partition, "A");
    assert_eq!(Some(1), node_a.observed);
    assert_abs_diff_eq!(Array1::from_vec(vec![0.0, 1.0]), node_a.profile, epsilon = 1e-12);

    let node_b = helpers::get_node_data(&graph, &partition, "B");
    assert_eq!(None, node_b.observed);
    assert_abs_diff_eq!(Array1::from_vec(vec![0.5, 0.5]), node_b.profile, epsilon = 1e-12);

    assert_eq!(graph.get_edges().len(), partition.edges.len());

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_rejects_tree_leaf_missing_from_metadata() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = BTreeMap::from([("A".to_owned(), "usa".to_owned())]);

    let error = attach_traits(&mut partition, &graph, &traits).unwrap_err().to_string();

    assert!(error.contains("tree leaves missing from metadata"), "{error}");
    assert!(error.contains('B'), "{error}");

    Ok(())
  }

  #[test]
  fn test_discrete_marginal_attach_traits_rejects_metadata_name_missing_from_tree() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let mut partition = helpers::make_partition(["usa", "germany"])?;
    let traits = BTreeMap::from([
      ("A".to_owned(), "usa".to_owned()),
      ("B".to_owned(), "germany".to_owned()),
      ("C".to_owned(), "usa".to_owned()),
    ]);

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
