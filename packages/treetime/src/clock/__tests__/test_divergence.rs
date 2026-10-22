#[cfg(test)]
mod tests {
  use crate::clock::divergence::{root_to_node_divergences, root_to_node_divergences_where_known};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::schedule::GraphSchedule;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::pretty_assert_abs_diff_eq;
  use treetime_utils::{assert_error, pretty_assert_map_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::all_named(          "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;", &[("root", 0.0), ("AB", 0.1), ("CD", 0.05), ("A", 0.2), ("B", 0.3), ("C", 0.25), ("D", 0.17)])]
  #[case::single_node(        "A:0.5;",                                                 &[("A", 0.0)])]
  #[case::linear_chain(       "((A:0.1)B:0.2)C:0.3;",                                   &[("C", 0.0), ("B", 0.2), ("A", 0.3)])]
  #[case::zero_branch_lengths("((A:0.0,B:0.1)AB:0.0,(C:0.2,D:0.0)CD:0.1)root:0.0;",   &[("root", 0.0), ("AB", 0.0), ("CD", 0.1), ("A", 0.0), ("B", 0.1), ("C", 0.3), ("D", 0.1)])]
  #[case::missing_lengths(    "((A,B:0.1)AB:0.2,C)root;",                               &[("root", 0.0), ("AB", 0.2), ("A", 0.2), ("B", 0.3), ("C", 0.0)])]
  #[trace]
  fn test_root_to_node_divergences_sum_branch_lengths_from_the_root(
    #[case] nwk: &str,
    #[case] expected: &[(&str, f64)],
  ) -> Result<(), Report> {
    let parsed = nwk_read(nwk.as_bytes())?;
    let names = parsed.names.clone();

    let actual: BTreeMap<String, f64> = root_to_node_divergences(&parsed.graph, |edge_key| parsed.branch_lengths[&edge_key].unwrap_or_default())?
      .into_iter()
      .map(|(key, div)| (names[key].to_owned(), div))
      .collect();

    let expected: BTreeMap<String, f64> = expected.iter().map(|(name, div)| ((*name).to_owned(), *div)).collect();
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::all_known(         "((A:0.125,B:0.25)AB:0.5,C:0.25)root;", &[("root", Some(0.0)), ("AB", Some(0.5)), ("A", Some(0.625)), ("B", Some(0.75)), ("C", Some(0.25))])]
  #[case::missing_leaf_edge( "((A,B:0.25)AB:0.5,C:0.25)root;",       &[("root", Some(0.0)), ("AB", Some(0.5)), ("A", None),        ("B", Some(0.75)), ("C", Some(0.25))])]
  #[case::missing_inner_edge("((A:0.125,B:0.25)AB,C:0.25)root;",     &[("root", Some(0.0)), ("AB", None),      ("A", None),        ("B", None),       ("C", Some(0.25))])]
  #[case::single_node(       "A;",                                   &[("A", Some(0.0))])]
  #[trace]
  fn test_root_to_node_divergences_where_known_leaves_nodes_below_a_missing_length_unknown(
    #[case] nwk: &str,
    #[case] expected: &[(&str, Option<f64>)],
  ) -> Result<(), Report> {
    let parsed = nwk_read(nwk.as_bytes())?;
    let names = parsed.names.clone();

    let actual: BTreeMap<String, Option<f64>> =
      root_to_node_divergences_where_known(&GraphSchedule::new(&parsed.graph)?, &parsed.branch_lengths)?
        .into_iter()
        .map(|(key, div)| (names[key].to_owned(), div))
        .collect();

    let expected: BTreeMap<String, Option<f64>> = expected.iter().map(|(name, div)| ((*name).to_owned(), *div)).collect();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_root_to_node_divergences_where_known_rejects_node_with_two_parents() -> Result<(), Report> {
    let mut graph = Graph::new();
    let root = graph.add_node();
    let left = graph.add_node();
    let right = graph.add_node();
    let shared = graph.add_node();
    let edges = [
      graph.add_edge(root, left)?,
      graph.add_edge(root, right)?,
      graph.add_edge(left, shared)?,
      graph.add_edge(right, shared)?,
    ];
    graph.build()?;
    let branch_lengths = edges.into_iter().map(|edge| (edge, Some(0.1))).collect();

    let result = root_to_node_divergences_where_known(&GraphSchedule::new(&graph)?, &branch_lengths);

    assert_error!(
      result,
      format!("The graph is not a tree: node {shared} has 2 parents, but a node of a tree has at most one")
    );
    Ok(())
  }

  #[test]
  fn test_root_to_node_divergences_gives_internal_nodes_their_path_sum() -> Result<(), Report> {
    let parsed = nwk_read(b"((A:0.1,B:0.2):0.3,C:0.4);".as_slice())?;
    let names = parsed.names.clone();

    let actual = root_to_node_divergences(&parsed.graph, |edge_key| {
      parsed.branch_lengths[&edge_key].unwrap_or_default()
    })?;

    let expected: BTreeMap<_, f64> = parsed
      .graph
      .get_nodes()
      .map(|node| {
        let div = match (node.is_root(), node.is_leaf(), &names[node.key()]) {
          (true, _, _) => 0.0,
          (false, false, _) => 0.3,
          (false, true, "A" | "C") => 0.4,
          (false, true, "B") => 0.5,
          (false, true, other) => panic!("unexpected leaf {other:?}"),
        };
        (node.key(), div)
      })
      .collect();
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_root_to_node_divergences_deep_tree() -> Result<(), Report> {
    let nwk = format!(
      "{};",
      (1..20).fold("A:0.05".to_owned(), |nwk, _| format!("({nwk}):0.05"))
    );
    let parsed = nwk_read(nwk.as_bytes())?;
    let names = parsed.names.clone();

    let actual = root_to_node_divergences(&parsed.graph, |edge_key| {
      parsed.branch_lengths[&edge_key].unwrap_or_default()
    })?;

    let leaf = parsed.graph.get_leaves().next().unwrap().key();
    assert_eq!("A", &names[leaf]);
    pretty_assert_abs_diff_eq!(19.0 * 0.05, actual[&leaf], epsilon = 1e-12);
    Ok(())
  }
}
