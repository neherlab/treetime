#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use eyre::Report;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.1)AB:0.1,C:0.1)root;";

  #[test]
  fn test_bad_leaves_flags_clock_outliers_among_dated_leaves() -> Result<(), Report> {
    let (graph, names) = helpers::tree()?;
    let constraints = helpers::dated(&graph, &names, &["B", "C"]);
    let outlier = find_node_key_by_name(&graph, &names, "B").expect("fixture node must exist");

    let actual = helpers::by_name(&names, &bad_leaves(&graph, &constraints, &btreeset! { outlier }));

    let expected = btreemap! {
      "A".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => false,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_bad_branches_undated_leaves_marks_leaves_without_a_date() -> Result<(), Report> {
    let (graph, names) = helpers::tree()?;
    let constraints = helpers::dated(&graph, &names, &["C", "AB"]);

    let actual = helpers::by_name(&names, &bad_leaves(&graph, &constraints, &BTreeSet::new()));

    let expected = btreemap! {
      "A".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => false,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_bad_branches_internal_node_is_bad_when_all_children_are_bad() -> Result<(), Report> {
    let (graph, names) = helpers::tree()?;
    let constraints = helpers::dated(&graph, &names, &["C"]);

    let bad_branches = derive_bad_branches(
      &graph,
      &constraints,
      &bad_leaves(&graph, &constraints, &BTreeSet::new()),
    )?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, helpers::by_name(&names, &bad_branches));
    Ok(())
  }

  #[test]
  fn test_bad_branches_internal_node_with_its_own_date_is_never_bad() -> Result<(), Report> {
    let (graph, names) = helpers::tree()?;
    let constraints = helpers::dated(&graph, &names, &["C", "AB"]);

    let bad_branches = derive_bad_branches(
      &graph,
      &constraints,
      &bad_leaves(&graph, &constraints, &BTreeSet::new()),
    )?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => false,
      "B".to_owned() => true,
      "C".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, helpers::by_name(&names, &bad_branches));
    Ok(())
  }

  #[test]
  fn test_bad_branches_dated_leaf_flagged_as_outlier_stays_bad() -> Result<(), Report> {
    let (graph, names) = helpers::tree()?;
    let constraints = helpers::dated(&graph, &names, &["A", "B", "C"]);
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
    let leaf_bad_branches = btreemap! {
      key("A") => true,
      key("B") => true,
      key("C") => false,
    };

    let bad_branches = derive_bad_branches(&graph, &constraints, &leaf_bad_branches)?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, helpers::by_name(&names, &bad_branches));
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn tree() -> Result<(Graph, BTreeMap<GraphNodeKey, Option<String>>), Report> {
      let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
      let names = nwk_parsed.names();
      Ok((nwk_parsed.graph, names))
    }

    pub(super) fn dated(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      dated: &[&str],
    ) -> DateConstraints {
      let date_constraints = dated
        .iter()
        .map(|name| {
          let key = find_node_key_by_name(graph, names, name).expect("fixture node must exist");
          (key, Some(Arc::new(Distribution::point(2020.0, 0.0))))
        })
        .collect();
      DateConstraints { date_constraints }
    }

    pub(super) fn by_name(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      flags: &BTreeMap<GraphNodeKey, bool>,
    ) -> BTreeMap<String, bool> {
      flags
        .iter()
        .map(|(key, bad)| (names[key].clone().expect("every fixture node is named"), *bad))
        .collect()
    }
  }
}
