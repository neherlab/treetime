#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::o;
  use crate::progress::NoopProgress;
  use crate::test_utils::dates_by_node;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::collections::BTreeSet;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::dates_csv::{DateConstraint, DateRange, DateValue};
  use treetime_io::nwk::nwk_read;

  #[test]
  fn test_load_date_constraints_success_three_leaves() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.15)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(2020.0),
      o!("B") => exact(2020.5),
      o!("C") => exact(2020.75),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::point(2020.0, 0.0))),
      loaded("B", Some(Distribution::point(2020.5, 0.0))),
      loaded("C", Some(Distribution::point(2020.75, 0.0))),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_mixed_leaves() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.15,D:0.18)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(2020.0),
      o!("B") => exact(2020.5),
      o!("C") => exact(2020.75),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::point(2020.0, 0.0))),
      loaded("B", Some(Distribution::point(2020.5, 0.0))),
      loaded("C", Some(Distribution::point(2020.75, 0.0))),
      loaded("D", None),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_range() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.15)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => range(2020.0, 2020.25),
      o!("B") => exact(2020.5),
      o!("C") => exact(2020.75),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::range((2020.0, 2020.25), 0.0))),
      loaded("B", Some(Distribution::point(2020.5, 0.0))),
      loaded("C", Some(Distribution::point(2020.75, 0.0))),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_internal_node() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,C:0.15)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(2020.0),
      o!("B") => exact(2020.5),
      o!("C") => exact(2020.75),
      o!("AB") => exact(2019.5),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::point(2020.0, 0.0))),
      loaded("AB", Some(Distribution::point(2019.5, 0.0))),
      loaded("B", Some(Distribution::point(2020.5, 0.0))),
      loaded("C", Some(Distribution::point(2020.75, 0.0))),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_derive_bad_branches_over_loaded_dates_marks_only_the_undated_leaf() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,(C:0.15,D:0.18)CD:0.1)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(2020.0),
      o!("B") => exact(2020.5),
      o!("C") => exact(2020.75),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::point(2020.0, 0.0))),
      loaded("AB", None),
      loaded("B", Some(Distribution::point(2020.5, 0.0))),
      loaded("C", Some(Distribution::point(2020.75, 0.0))),
      loaded("CD", None),
      loaded("D", None),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    let expected_bad_branches = btreemap! {
      o!("A") => false,
      o!("AB") => false,
      o!("B") => false,
      o!("C") => false,
      o!("CD") => false,
      o!("D") => true,
      o!("root") => false,
    };
    assert_eq!(
      expected_bad_branches,
      derived_bad_branches(&names, &graph, &constraints)?
    );
    Ok(())
  }

  #[test]
  fn test_derive_bad_branches_over_loaded_dates_marks_a_subtree_whose_children_are_all_bad() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(((A:0.1,B:0.2)AB:0.1,(C:0.15,D:0.18)CD:0.1)ABCD:0.1,E:0.2)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("C") => exact(2020.0),
      o!("D") => exact(2020.5),
      o!("E") => exact(2020.75),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", None),
      loaded("AB", None),
      loaded("ABCD", None),
      loaded("B", None),
      loaded("C", Some(Distribution::point(2020.0, 0.0))),
      loaded("CD", None),
      loaded("D", Some(Distribution::point(2020.5, 0.0))),
      loaded("E", Some(Distribution::point(2020.75, 0.0))),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    let expected_bad_branches = btreemap! {
      o!("A") => true,
      o!("AB") => true,
      o!("ABCD") => false,
      o!("B") => true,
      o!("C") => false,
      o!("CD") => false,
      o!("D") => false,
      o!("E") => false,
      o!("root") => false,
    };
    assert_eq!(
      expected_bad_branches,
      derived_bad_branches(&names, &graph, &constraints)?
    );
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_boundary_exactly_three_leaves() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.15)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(2020.0),
      o!("B") => exact(2020.5),
      o!("C") => exact(2020.75),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    assert_eq!(actual.iter().filter(|n| n.date_constraint.is_some()).count(), 3);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_date_with_none_value() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.15,D:0.18)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(2020.0),
      o!("B") => None,
      o!("C") => exact(2020.5),
      o!("D") => exact(2020.75),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::point(2020.0, 0.0))),
      loaded("B", None),
      loaded("C", Some(Distribution::point(2020.5, 0.0))),
      loaded("D", Some(Distribution::point(2020.75, 0.0))),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_derive_bad_branches_over_loaded_dates_stops_at_the_first_dated_descendant() -> Result<(), Report> {
    let nwk_parsed = nwk_read(
      b"((((((A:0.1,B:0.1)L1:0.1,C:0.1)L2:0.1,D:0.1)L3:0.1,E:0.1)L4:0.1,F:0.1)L5:0.1,G:0.1)root:0.0;".as_slice(),
    )?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("C") => exact(2020.0),
      o!("D") => exact(2020.25),
      o!("E") => exact(2020.5),
      o!("F") => exact(2020.75),
      o!("G") => exact(2021.0),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", None),
      loaded("B", None),
      loaded("C", Some(Distribution::point(2020.0, 0.0))),
      loaded("D", Some(Distribution::point(2020.25, 0.0))),
      loaded("E", Some(Distribution::point(2020.5, 0.0))),
      loaded("F", Some(Distribution::point(2020.75, 0.0))),
      loaded("G", Some(Distribution::point(2021.0, 0.0))),
      loaded("L1", None),
      loaded("L2", None),
      loaded("L3", None),
      loaded("L4", None),
      loaded("L5", None),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    let expected_bad_branches = btreemap! {
      o!("A") => true,
      o!("B") => true,
      o!("C") => false,
      o!("D") => false,
      o!("E") => false,
      o!("F") => false,
      o!("G") => false,
      o!("L1") => true,
      o!("L2") => false,
      o!("L3") => false,
      o!("L4") => false,
      o!("L5") => false,
      o!("root") => false,
    };
    assert_eq!(
      expected_bad_branches,
      derived_bad_branches(&names, &graph, &constraints)?
    );
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_wide_tree() -> Result<(), Report> {
    let nwk_parsed =
      nwk_read(b"(A:0.1,B:0.1,C:0.1,D:0.1,E:0.1,F:0.1,G:0.1,H:0.1,I:0.1,J:0.1,K:0.1,L:0.1)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(2020.0),
      o!("C") => exact(2020.2),
      o!("E") => exact(2020.4),
      o!("G") => exact(2020.6),
      o!("I") => exact(2020.8),
      o!("K") => exact(2021.0),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::point(2020.0, 0.0))),
      loaded("B", None),
      loaded("C", Some(Distribution::point(2020.2, 0.0))),
      loaded("D", None),
      loaded("E", Some(Distribution::point(2020.4, 0.0))),
      loaded("F", None),
      loaded("G", Some(Distribution::point(2020.6, 0.0))),
      loaded("H", None),
      loaded("I", Some(Distribution::point(2020.8, 0.0))),
      loaded("J", None),
      loaded("K", Some(Distribution::point(2021.0, 0.0))),
      loaded("L", None),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_mixed_ranges_and_points() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.15,D:0.18)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => range(2020.0, 2020.25),
      o!("B") => exact(2020.5),
      o!("C") => range(2020.6, 2020.8),
      o!("D") => exact(2021.0),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::range((2020.0, 2020.25), 0.0))),
      loaded("B", Some(Distribution::point(2020.5, 0.0))),
      loaded("C", Some(Distribution::range((2020.6, 2020.8), 0.0))),
      loaded("D", Some(Distribution::point(2021.0, 0.0))),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_internal_node_with_range() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,C:0.15)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(2020.0),
      o!("B") => exact(2020.5),
      o!("C") => exact(2020.75),
      o!("AB") => range(2019.0, 2019.75),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::point(2020.0, 0.0))),
      loaded("AB", Some(Distribution::range((2019.0, 2019.75), 0.0))),
      loaded("B", Some(Distribution::point(2020.5, 0.0))),
      loaded("C", Some(Distribution::point(2020.75, 0.0))),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_negative_time() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"(A:0.1,B:0.2,C:0.15)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let dates: BTreeMap<String, Option<DateConstraint>> = btreemap! {
      o!("A") => exact(-500.0),
      o!("B") => exact(-250.0),
      o!("C") => exact(0.0),
    };

    let constraints = load_date_constraints(&dates_by_node(dates, &graph, &names), &graph, &NoopProgress)?;

    let actual = node_constraints(&names, &graph, &constraints);
    let expected = vec![
      loaded("A", Some(Distribution::point(-500.0, 0.0))),
      loaded("B", Some(Distribution::point(-250.0, 0.0))),
      loaded("C", Some(Distribution::point(0.0, 0.0))),
      loaded("root", None),
    ];
    assert_eq!(actual, expected);
    Ok(())
  }

  mod helpers {
    use super::*;

    #[derive(Clone, Default, Debug, PartialEq)]
    pub(super) struct LoadedNode {
      pub(super) name: Option<String>,
      pub(super) date_constraint: Option<Arc<Distribution<NegLog>>>,
    }

    pub(super) fn node_constraints(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      constraints: &DateConstraints,
    ) -> Vec<LoadedNode> {
      graph
        .get_nodes()
        .map(|node| LoadedNode {
          name: names.get(&node.key()).cloned().flatten(),
          date_constraint: constraints.by_node[&node.key()].clone(),
        })
        .sorted_by_key(|n| n.name.clone().unwrap_or_default())
        .collect_vec()
    }

    pub(super) fn derived_bad_branches(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      constraints: &DateConstraints,
    ) -> Result<BTreeMap<String, bool>, Report> {
      let leaf_bad_branches = bad_leaves(graph, constraints, &BTreeSet::new());
      Ok(
        derive_bad_branches(graph, constraints, &leaf_bad_branches)?
          .into_iter()
          .map(|(key, bad)| (names[&key].clone().expect("every fixture node is named"), bad))
          .collect(),
      )
    }

    pub(super) fn loaded(name: &str, date_constraint: Option<Distribution<NegLog>>) -> LoadedNode {
      LoadedNode {
        name: Some(name.to_owned()),
        date_constraint: date_constraint.map(Arc::new),
      }
    }

    pub(super) fn exact(value: f64) -> Option<DateConstraint> {
      Some(DateConstraint::exact(value))
    }

    pub(super) fn range(start: f64, end: f64) -> Option<DateConstraint> {
      Some(DateConstraint {
        raw: format!("{start}/{end}"),
        value: DateValue::Range(DateRange { start, end }),
      })
    }
  }

  use helpers::*;
}
