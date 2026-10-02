#[cfg(test)]
mod tests {
  use crate::clock::divergence::root_to_node_divergences;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

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
    let parsed = nwk_read_str(nwk)?;
    let names = parsed.names();

    let actual: BTreeMap<String, f64> = root_to_node_divergences(&parsed.graph, |edge_key| parsed.branch_lengths[&edge_key].unwrap_or_default())?
      .into_iter()
      .map(|(key, div)| (names[&key].clone().unwrap(), div))
      .collect();

    let expected: BTreeMap<String, f64> = expected.iter().map(|(name, div)| ((*name).to_owned(), *div)).collect();
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_root_to_node_divergences_covers_unnamed_nodes() -> Result<(), Report> {
    let parsed = nwk_read_str("((A:0.1,B:0.2):0.1,(C:0.2,D:0.12):0.05):0.01;")?;

    let actual = root_to_node_divergences(&parsed.graph, |edge_key| {
      parsed.branch_lengths[&edge_key].unwrap_or_default()
    })?;

    let keys = parsed.graph.get_nodes().map(|node| node.key()).collect::<Vec<_>>();
    assert_eq!(keys, actual.keys().copied().collect::<Vec<_>>());
    Ok(())
  }

  #[test]
  fn test_root_to_node_divergences_deep_tree() -> Result<(), Report> {
    let nwk = format!(
      "{};",
      (1..20).fold("A:0.05".to_owned(), |nwk, _| format!("({nwk}):0.05"))
    );
    let parsed = nwk_read_str(&nwk)?;
    let names = parsed.names();

    let actual = root_to_node_divergences(&parsed.graph, |edge_key| {
      parsed.branch_lengths[&edge_key].unwrap_or_default()
    })?;

    let leaf = parsed.graph.get_leaves().next().unwrap().key();
    assert_eq!(Some("A"), names[&leaf].as_deref());
    assert_abs_diff_eq!(19.0 * 0.05, actual[&leaf], epsilon = 1e-12);
    Ok(())
  }
}
