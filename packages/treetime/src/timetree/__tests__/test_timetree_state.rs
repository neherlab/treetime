#[cfg(test)]
mod tests {
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use crate::clock::date_constraints::DateConstraints;
  use crate::test_utils::find_node_key_by_name;
  use ndarray::array;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_utils::assert_error;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_timetree_state_reset_date_edges_clears_distribution_and_message() -> Result<(), Report> {
    let graph = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?.graph;
    let mut state = TimetreeState::new(&graph);
    for edge_ref in graph.get_edges() {
      let key = edge_ref.key();
      let entry = state.edge_mut(key);
      entry.branch_length_distribution = Some(Arc::new(Distribution::point(1.0, 0.0)));
      entry.msg_to_parent = Some(Arc::new(Distribution::point(2.0, 0.0)));
    }

    state.reset_date_edges_for_topology_change(&graph);

    for edge_ref in graph.get_edges() {
      let key = edge_ref.key();
      let entry = state.edge(key);
      assert_eq!(None, entry.branch_length_distribution);
      assert_eq!(None, entry.msg_to_parent);
    }

    Ok(())
  }

  #[test]
  fn test_timetree_state_reseed_preserves_distribution_and_message() -> Result<(), Report> {
    let graph = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?.graph;
    let mut state = TimetreeState::new(&graph);
    let key = graph
      .get_edges()
      .collect::<Vec<_>>()
      .first()
      .expect("tree has at least one edge")
      .key();
    let dist = Arc::new(Distribution::point(3.0, 0.0));
    let msg = Arc::new(Distribution::point(4.0, 0.0));
    {
      let entry = state.edge_mut(key);
      entry.branch_length_distribution = Some(Arc::clone(&dist));
      entry.msg_to_parent = Some(Arc::clone(&msg));
    }

    state.reseed_from_values(&graph);

    let entry = state.edge(key);
    assert_eq!(Some(dist), entry.branch_length_distribution);
    assert_eq!(Some(msg), entry.msg_to_parent);

    Ok(())
  }

  #[test]
  fn test_timetree_state_likely_times_names_node_with_nan_distribution() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:1.0,B:1.0)I:1.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = find_node_key_by_name(&graph, &names, "I").expect("internal node I not found");
    let mut state = TimetreeState::new(&graph);
    state.node_mut(key).time_distribution = Some(Arc::new(helpers::nan_distribution()?));

    let expected = format!(
      "When finding the most likely time of node {key}: \
       Cannot find the most likely time of a distribution function: its values contain NaN"
    );
    assert_error!(state.likely_times(&DateConstraints::default()), expected);
    assert_error!(state.coalescent_node_times(), expected);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn nan_distribution() -> Result<Distribution<NegLog>, Report> {
      Distribution::function(array![0.0, 1.0, 2.0], array![1.0, f64::NAN, 2.0])
    }
  }
}
