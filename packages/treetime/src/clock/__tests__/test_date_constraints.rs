#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::o;
  use crate::progress::NoopProgress;
  use crate::test_utils::dates_by_node;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use deser::{Deserialize, Serialize};
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
  use treetime_utils::io::json::json_read_str;

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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "B", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "C", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "B", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "C", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "D", "date_constraint": null},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"range": {"range": [2020.0, 2020.25], "ampl": 0.0}}},
        {"name": "B", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "C", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "AB", "date_constraint": {"point": {"t": 2019.5, "ampl": 0.0}}},
        {"name": "B", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "C", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "AB", "date_constraint": null},
        {"name": "B", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "C", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "CD", "date_constraint": null},
        {"name": "D", "date_constraint": null},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": null},
        {"name": "AB", "date_constraint": null},
        {"name": "ABCD", "date_constraint": null},
        {"name": "B", "date_constraint": null},
        {"name": "C", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "CD", "date_constraint": null},
        {"name": "D", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "E", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "B", "date_constraint": null},
        {"name": "C", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "D", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": null},
        {"name": "B", "date_constraint": null},
        {"name": "C", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "D", "date_constraint": {"point": {"t": 2020.25, "ampl": 0.0}}},
        {"name": "E", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "F", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "G", "date_constraint": {"point": {"t": 2021.0, "ampl": 0.0}}},
        {"name": "L1", "date_constraint": null},
        {"name": "L2", "date_constraint": null},
        {"name": "L3", "date_constraint": null},
        {"name": "L4", "date_constraint": null},
        {"name": "L5", "date_constraint": null},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "B", "date_constraint": null},
        {"name": "C", "date_constraint": {"point": {"t": 2020.2, "ampl": 0.0}}},
        {"name": "D", "date_constraint": null},
        {"name": "E", "date_constraint": {"point": {"t": 2020.4, "ampl": 0.0}}},
        {"name": "F", "date_constraint": null},
        {"name": "G", "date_constraint": {"point": {"t": 2020.6, "ampl": 0.0}}},
        {"name": "H", "date_constraint": null},
        {"name": "I", "date_constraint": {"point": {"t": 2020.8, "ampl": 0.0}}},
        {"name": "J", "date_constraint": null},
        {"name": "K", "date_constraint": {"point": {"t": 2021.0, "ampl": 0.0}}},
        {"name": "L", "date_constraint": null},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"range": {"range": [2020.0, 2020.25], "ampl": 0.0}}},
        {"name": "B", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "C", "date_constraint": {"range": {"range": [2020.6, 2020.8], "ampl": 0.0}}},
        {"name": "D", "date_constraint": {"point": {"t": 2021.0, "ampl": 0.0}}},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"point": {"t": 2020.0, "ampl": 0.0}}},
        {"name": "AB", "date_constraint": {"range": {"range": [2019.0, 2019.75], "ampl": 0.0}}},
        {"name": "B", "date_constraint": {"point": {"t": 2020.5, "ampl": 0.0}}},
        {"name": "C", "date_constraint": {"point": {"t": 2020.75, "ampl": 0.0}}},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
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
    let expected: Vec<LoadedNode> = json_read_str(
      r#"[
        {"name": "A", "date_constraint": {"point": {"t": -500.0, "ampl": 0.0}}},
        {"name": "B", "date_constraint": {"point": {"t": -250.0, "ampl": 0.0}}},
        {"name": "C", "date_constraint": {"point": {"t": 0.0, "ampl": 0.0}}},
        {"name": "root", "date_constraint": null}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  mod helpers {
    use super::*;

    #[derive(Clone, Default, Debug, PartialEq, Serialize, Deserialize)]
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
