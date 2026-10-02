#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::test_utils::{find_node_key_by_name, point_date_constraints};
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use eyre::Report;
  use generators::gen_bad_branch_case;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use std::collections::{BTreeMap, BTreeSet};
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.1)AB:0.1,C:0.1)root;";

  const DATE: f64 = 2020.0;

  #[test]
  fn test_bad_leaves_flags_clock_outliers_among_dated_leaves() -> Result<(), Report> {
    let (graph, names) = helpers::tree(TREE_NEWICK)?;
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
    let (graph, names) = helpers::tree(TREE_NEWICK)?;
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
    let actual = helpers::derived_by_name(TREE_NEWICK, &["C"], &[])?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_bad_branches_internal_node_with_its_own_date_is_never_bad() -> Result<(), Report> {
    let actual = helpers::derived_by_name(TREE_NEWICK, &["C", "AB"], &[])?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => false,
      "B".to_owned() => true,
      "C".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_bad_branches_dated_leaf_flagged_as_outlier_stays_bad() -> Result<(), Report> {
    let (graph, names) = helpers::tree(TREE_NEWICK)?;
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

  #[test]
  fn test_bad_branches_multifurcation_is_good_when_one_of_its_children_is_good() -> Result<(), Report> {
    let actual = helpers::derived_by_name("((A:0.1,B:0.1,C:0.1)P:0.1,D:0.1)root;", &["C", "D"], &[])?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => false,
      "D".to_owned() => false,
      "P".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_bad_branches_multifurcation_is_bad_when_all_of_its_children_are_bad() -> Result<(), Report> {
    let actual = helpers::derived_by_name("((A:0.1,B:0.1,C:0.1)P:0.1,D:0.1)root;", &["D"], &[])?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => true,
      "D".to_owned() => false,
      "P".to_owned() => true,
      "root".to_owned() => false,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_bad_branches_root_is_bad_when_nothing_below_it_is_dated() -> Result<(), Report> {
    let actual = helpers::derived_by_name(TREE_NEWICK, &[], &[])?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => true,
      "root".to_owned() => true,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_bad_branches_nested_dated_nodes_keep_every_ancestor_good() -> Result<(), Report> {
    let actual = helpers::derived_by_name("(((A:0.1,B:0.1)X:0.1,C:0.1)Y:0.1,D:0.1)root;", &["X", "Y"], &[])?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => true,
      "D".to_owned() => true,
      "X".to_owned() => false,
      "Y".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_bad_branches_internal_node_with_a_date_range_is_never_bad() -> Result<(), Report> {
    let actual = helpers::derived_by_name(TREE_NEWICK, &[], &["AB"])?;

    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => false,
      "B".to_owned() => true,
      "C".to_owned() => true,
      "root".to_owned() => false,
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_prop_bad_branches_node_is_bad_iff_its_subtree_holds_no_evidence(case in gen_bad_branch_case()) {
      let (graph, names) = helpers::tree(&case.newick).unwrap();
      let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("generated node must exist");
      let constraints = DateConstraints {
        by_node: case
          .dated_internal
          .iter()
          .map(|name| (key(name), Some(Arc::new(Distribution::point(DATE, 0.0)))))
          .collect(),
      };
      let leaf_bad_branches: BTreeMap<GraphNodeKey, bool> =
        case.bad_leaves.iter().map(|(name, bad)| (key(name), *bad)).collect();
      let has_evidence = |node_key: GraphNodeKey| {
        let node = graph.get_node(node_key).expect("node must exist");
        if node.is_leaf() {
          !leaf_bad_branches[&node_key]
        } else {
          constraints.date_constraint(node_key).is_some()
        }
      };

      let actual = derive_bad_branches(&graph, &constraints, &leaf_bad_branches).unwrap();

      let expected: BTreeMap<GraphNodeKey, bool> = graph
        .get_nodes()
        .map(|node| {
          let key = node.key();
          (key, !helpers::subtree_keys(&graph, key).into_iter().any(has_evidence))
        })
        .collect();
      prop_assert_eq!(expected, actual);
    }
  }

  mod generators {
    use super::*;

    #[derive(Debug, Clone)]
    pub(super) struct BadBranchCase {
      pub newick: String,
      pub bad_leaves: Vec<(String, bool)>,
      pub dated_internal: Vec<String>,
    }

    pub(super) fn gen_bad_branch_case() -> impl Strategy<Value = BadBranchCase> {
      (2_usize..=9).prop_flat_map(|n_leaves| {
        (
          prop::collection::vec(
            (
              any::<prop::sample::Index>(),
              any::<prop::sample::Index>(),
              any::<prop::sample::Index>(),
              any::<bool>(),
            ),
            n_leaves - 1,
          ),
          prop::collection::vec(any::<bool>(), n_leaves),
          prop::collection::vec(any::<bool>(), n_leaves - 1),
        )
          .prop_map(move |(merges, leaf_flags, internal_flags)| {
            let leaves = (0..n_leaves).map(|index| format!("L{index}")).collect::<Vec<_>>();
            let (newick, internal) = random_newick(&leaves, &merges);
            BadBranchCase {
              newick,
              bad_leaves: leaves.into_iter().zip(leaf_flags).collect(),
              dated_internal: internal
                .into_iter()
                .zip(internal_flags)
                .filter_map(|(name, dated)| dated.then_some(name))
                .collect(),
            }
          })
      })
    }

    fn random_newick(
      leaves: &[String],
      merges: &[(prop::sample::Index, prop::sample::Index, prop::sample::Index, bool)],
    ) -> (String, Vec<String>) {
      let mut subtrees = leaves.to_vec();
      let mut internal = vec![];
      for (index, (first, second, third, multifurcate)) in merges.iter().enumerate() {
        if subtrees.len() < 2 {
          break;
        }
        let mut children = vec![
          subtrees.remove(first.index(subtrees.len())),
          subtrees.remove(second.index(subtrees.len())),
        ];
        if *multifurcate && !subtrees.is_empty() {
          children.push(subtrees.remove(third.index(subtrees.len())));
        }
        let name = format!("I{index}");
        let children = children.iter().map(|child| format!("{child}:0.1")).collect::<Vec<_>>();
        subtrees.push(format!("({}){name}", children.join(",")));
        internal.push(name);
      }
      (format!("{};", subtrees.concat()), internal)
    }
  }

  mod helpers {
    use super::*;

    pub(super) fn tree(newick: &str) -> Result<(Graph, BTreeMap<GraphNodeKey, Option<String>>), Report> {
      let nwk_parsed = nwk_read_str(newick)?;
      let names = nwk_parsed.names();
      Ok((nwk_parsed.graph, names))
    }

    pub(super) fn dated(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      dated: &[&str],
    ) -> DateConstraints {
      let dates = dated.iter().map(|name| (*name, DATE)).collect::<Vec<_>>();
      point_date_constraints(graph, names, &dates)
    }

    pub(super) fn derived_by_name(
      newick: &str,
      dated_names: &[&str],
      ranged_names: &[&str],
    ) -> Result<BTreeMap<String, bool>, Report> {
      let (graph, names) = tree(newick)?;
      let mut constraints = dated(&graph, &names, dated_names);
      for name in ranged_names {
        let key = find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");
        let range = Distribution::range((DATE - 1.0, DATE), 0.0);
        constraints.by_node.insert(key, Some(Arc::new(range)));
      }
      let bad_branches = derive_bad_branches(
        &graph,
        &constraints,
        &bad_leaves(&graph, &constraints, &BTreeSet::new()),
      )?;
      Ok(by_name(&names, &bad_branches))
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

    pub(super) fn subtree_keys(graph: &Graph, key: GraphNodeKey) -> Vec<GraphNodeKey> {
      let mut keys = vec![];
      let mut stack = vec![key];
      while let Some(current) = stack.pop() {
        keys.push(current);
        let node = graph.get_node(current).expect("node must exist");
        for edge_key in node.outbound() {
          stack.push(graph.get_edge(*edge_key).expect("edge must exist").target());
        }
      }
      keys
    }
  }
}
