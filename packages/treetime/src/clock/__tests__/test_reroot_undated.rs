#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_length_or_zero;
  use crate::clock::divergence::root_to_node_divergences;
  use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use crate::clock::reroot::reroot_undated;
  use crate::progress::NoopProgress;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read;

  fn root_to_tip_distances(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    let divergences = root_to_node_divergences(graph, |key| branch_length_or_zero(branch_lengths, key)).unwrap();
    let mut distances: Vec<f64> = graph.get_leaves().map(|leaf| divergences[&leaf.key()]).collect();
    distances.sort_by(f64::total_cmp);
    distances
  }

  #[test]
  fn test_reroot_undated_min_dev_roots_a_cherry_at_equal_distances() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.3)root;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    reroot_undated(
      &mut graph,
      &mut branch_lengths,
      &RerootSpec::Method(RerootMethod::MinDev),
      &names,
      &NoopProgress,
    )?;

    let dists = root_to_tip_distances(&graph, &branch_lengths);
    assert_eq!(2, dists.len());
    assert_abs_diff_eq!(dists[0], 0.2, epsilon = 1e-15);
    assert_abs_diff_eq!(dists[1], 0.2, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_reroot_undated_min_dev_places_the_root_at_the_analytic_minimum_of_the_variance() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:1,B:2,C:6)R;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;

    reroot_undated(
      &mut graph,
      &mut branch_lengths,
      &RerootSpec::Method(RerootMethod::MinDev),
      &names,
      &NoopProgress,
    )?;

    let dists = root_to_tip_distances(&graph, &branch_lengths);
    assert_abs_diff_eq!(dists.as_slice(), [3.25, 3.75, 4.25].as_slice(), epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_reroot_undated_min_dev_removes_an_undated_stem_before_the_root_search() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(((A:0.1,B:0.3)X:0.05,C:0.4)R:0.7)STEM;".as_slice())?;
    let names = nwk_parsed.names();
    let mut graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let unstemmed = nwk_read(b"((A:0.1,B:0.3)X:0.05,C:0.4)R;".as_slice())?;
    let unstemmed_names = unstemmed.names();
    let mut unstemmed_graph = unstemmed.graph;
    let mut unstemmed_branch_lengths = unstemmed.branch_lengths;

    let result = reroot_undated(
      &mut graph,
      &mut branch_lengths,
      &RerootSpec::Method(RerootMethod::MinDev),
      &names,
      &NoopProgress,
    )?;
    reroot_undated(
      &mut unstemmed_graph,
      &mut unstemmed_branch_lengths,
      &RerootSpec::Method(RerootMethod::MinDev),
      &unstemmed_names,
      &NoopProgress,
    )?;

    assert!(result.stem_removal.is_some());
    assert_eq!(3, graph.get_leaves().count());
    let expected = root_to_tip_distances(&unstemmed_graph, &unstemmed_branch_lengths);
    let actual = root_to_tip_distances(&graph, &branch_lengths);
    assert_abs_diff_eq!(expected.as_slice(), actual.as_slice(), epsilon = 1e-12);
    Ok(())
  }
}
