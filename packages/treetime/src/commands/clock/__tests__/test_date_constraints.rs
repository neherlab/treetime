#[cfg(test)]
mod tests {
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::o;
  use crate::representation::payload::traits::DateConstraintNode;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use serde::{Deserialize, Serialize};
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::edge::GraphEdge;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::{GraphNode, Named, TimeConstraint};
  use treetime_io::dates_csv::DatesMap;
  use treetime_io::nwk::{EdgeFromNwk, NodeFromNwk, nwk_read_str};
  use treetime_utils::io::json::json_read_str;

  #[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq)]
  struct TestNode {
    name: Option<String>,
    time_distribution: Option<Arc<Distribution>>,
    bad_branch: bool,
  }

  impl GraphNode for TestNode {}

  impl TimeConstraint<Arc<Distribution>> for TestNode {
    fn time_distribution(&self) -> &Option<Arc<Distribution>> {
      &self.time_distribution
    }
    fn set_time_distribution(&mut self, dist: Option<Arc<Distribution>>) {
      self.time_distribution = dist;
    }
    fn bad_branch(&self) -> bool {
      self.bad_branch
    }
    fn set_bad_branch(&mut self, bad: bool) {
      self.bad_branch = bad;
    }
  }

  impl DateConstraintNode for TestNode {}

  impl Named for TestNode {
    fn name(&self) -> Option<impl AsRef<str>> {
      self.name.as_deref()
    }
    fn set_name(&mut self, name: Option<impl AsRef<str>>) {
      self.name = name.map(|n| o!(n.as_ref()));
    }
  }

  impl NodeFromNwk for TestNode {
    fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
      Ok(Self {
        name: name.map(|n| o!(n.as_ref())),
        ..Default::default()
      })
    }
  }

  #[derive(Clone, Default, Debug, Serialize, Deserialize)]
  struct TestEdge;

  impl GraphEdge for TestEdge {}

  impl EdgeFromNwk for TestEdge {
    fn from_nwk(_: Option<f64>) -> Result<Self, Report> {
      Ok(Self)
    }
  }

  type TestGraph = Graph<TestNode, TestEdge, ()>;

  fn get_node_payloads(graph: &TestGraph) -> Vec<TestNode> {
    graph
      .get_node_payloads()
      .map(|node| node.read_arc().clone())
      .sorted_by_key(|n| n.name().map(|n| o!(n.as_ref())).unwrap_or_default())
      .collect_vec()
  }

  #[test]
  fn test_load_date_constraints_success_three_leaves() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(A:0.1,B:0.2,C:0.15)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFraction": 2020.75}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_mixed_leaves() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(A:0.1,B:0.2,C:0.15,D:0.18)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFraction": 2020.75}
      }"#,
    )?;
    // D has no date

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "D", "time_distribution": null, "bad_branch": true},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_range() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(A:0.1,B:0.2,C:0.15)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFractionRange": [2020.0, 2020.25]},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFraction": 2020.75}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Range": {"range": [2020.0, 2020.25], "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_internal_node() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.15)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFraction": 2020.75},
        "AB": {"YearFraction": 2019.5}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "AB", "time_distribution": {"Point": {"t": 2019.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_bad_branch_propagation() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.15,D:0.18)CD:0.1)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFraction": 2020.75}
      }"#,
    )?;
    // D has no date

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "AB", "time_distribution": null, "bad_branch": false},
        {"name": "B", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "CD", "time_distribution": null, "bad_branch": false},
        {"name": "D", "time_distribution": null, "bad_branch": true},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_all_children_bad() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(((A:0.1,B:0.2)AB:0.1,(C:0.15,D:0.18)CD:0.1)ABCD:0.1,E:0.2)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "C": {"YearFraction": 2020.0},
        "D": {"YearFraction": 2020.5},
        "E": {"YearFraction": 2020.75}
      }"#,
    )?;
    // A and B have no dates

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": null, "bad_branch": true},
        {"name": "AB", "time_distribution": null, "bad_branch": true},
        {"name": "ABCD", "time_distribution": null, "bad_branch": false},
        {"name": "B", "time_distribution": null, "bad_branch": true},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "CD", "time_distribution": null, "bad_branch": false},
        {"name": "D", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "E", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_boundary_exactly_three_leaves() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(A:0.1,B:0.2,C:0.15)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFraction": 2020.75}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    assert_eq!(actual.iter().filter(|n| n.time_distribution.is_some()).count(), 3);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_date_with_none_value() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(A:0.1,B:0.2,C:0.15,D:0.18)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "B": null,
        "C": {"YearFraction": 2020.5},
        "D": {"YearFraction": 2020.75}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": null, "bad_branch": true},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "D", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_deep_tree_propagation() -> Result<(), Report> {
    let graph: TestGraph =
      nwk_read_str("((((((A:0.1,B:0.1)L1:0.1,C:0.1)L2:0.1,D:0.1)L3:0.1,E:0.1)L4:0.1,F:0.1)L5:0.1,G:0.1)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "C": {"YearFraction": 2020.0},
        "D": {"YearFraction": 2020.25},
        "E": {"YearFraction": 2020.5},
        "F": {"YearFraction": 2020.75},
        "G": {"YearFraction": 2021.0}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": null, "bad_branch": true},
        {"name": "B", "time_distribution": null, "bad_branch": true},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "D", "time_distribution": {"Point": {"t": 2020.25, "ampl": 1.0}}, "bad_branch": false},
        {"name": "E", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "F", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "G", "time_distribution": {"Point": {"t": 2021.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "L1", "time_distribution": null, "bad_branch": true},
        {"name": "L2", "time_distribution": null, "bad_branch": false},
        {"name": "L3", "time_distribution": null, "bad_branch": false},
        {"name": "L4", "time_distribution": null, "bad_branch": false},
        {"name": "L5", "time_distribution": null, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_wide_tree() -> Result<(), Report> {
    let graph: TestGraph =
      nwk_read_str("(A:0.1,B:0.1,C:0.1,D:0.1,E:0.1,F:0.1,G:0.1,H:0.1,I:0.1,J:0.1,K:0.1,L:0.1)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "C": {"YearFraction": 2020.2},
        "E": {"YearFraction": 2020.4},
        "G": {"YearFraction": 2020.6},
        "I": {"YearFraction": 2020.8},
        "K": {"YearFraction": 2021.0}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": null, "bad_branch": true},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.2, "ampl": 1.0}}, "bad_branch": false},
        {"name": "D", "time_distribution": null, "bad_branch": true},
        {"name": "E", "time_distribution": {"Point": {"t": 2020.4, "ampl": 1.0}}, "bad_branch": false},
        {"name": "F", "time_distribution": null, "bad_branch": true},
        {"name": "G", "time_distribution": {"Point": {"t": 2020.6, "ampl": 1.0}}, "bad_branch": false},
        {"name": "H", "time_distribution": null, "bad_branch": true},
        {"name": "I", "time_distribution": {"Point": {"t": 2020.8, "ampl": 1.0}}, "bad_branch": false},
        {"name": "J", "time_distribution": null, "bad_branch": true},
        {"name": "K", "time_distribution": {"Point": {"t": 2021.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "L", "time_distribution": null, "bad_branch": true},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_mixed_ranges_and_points() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(A:0.1,B:0.2,C:0.15,D:0.18)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFractionRange": [2020.0, 2020.25]},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFractionRange": [2020.6, 2020.8]},
        "D": {"YearFraction": 2021.0}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Range": {"range": [2020.0, 2020.25], "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "C", "time_distribution": {"Range": {"range": [2020.6, 2020.8], "ampl": 1.0}}, "bad_branch": false},
        {"name": "D", "time_distribution": {"Point": {"t": 2021.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_idempotency() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(A:0.1,B:0.2,C:0.15)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFraction": 2020.75}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;
    let first_run = get_node_payloads(&graph);

    load_date_constraints(&dates, &graph)?;
    let second_run = get_node_payloads(&graph);

    assert_eq!(first_run, second_run);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_internal_node_with_range() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.15)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": 2020.0},
        "B": {"YearFraction": 2020.5},
        "C": {"YearFraction": 2020.75},
        "AB": {"YearFractionRange": [2019.0, 2019.75]}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Point": {"t": 2020.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "AB", "time_distribution": {"Range": {"range": [2019.0, 2019.75], "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": {"Point": {"t": 2020.5, "ampl": 1.0}}, "bad_branch": false},
        {"name": "C", "time_distribution": {"Point": {"t": 2020.75, "ampl": 1.0}}, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }

  #[test]
  fn test_load_date_constraints_negative_time() -> Result<(), Report> {
    let graph: TestGraph = nwk_read_str("(A:0.1,B:0.2,C:0.15)root:0.0;")?;
    let dates: DatesMap = json_read_str(
      r#"{
        "A": {"YearFraction": -500.0},
        "B": {"YearFraction": -250.0},
        "C": {"YearFraction": 0.0}
      }"#,
    )?;

    load_date_constraints(&dates, &graph)?;

    let actual = get_node_payloads(&graph);
    let expected: Vec<TestNode> = json_read_str(
      r#"[
        {"name": "A", "time_distribution": {"Point": {"t": -500.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "B", "time_distribution": {"Point": {"t": -250.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "C", "time_distribution": {"Point": {"t": 0.0, "ampl": 1.0}}, "bad_branch": false},
        {"name": "root", "time_distribution": null, "bad_branch": false}
      ]"#,
    )?;
    assert_eq!(actual, expected);
    Ok(())
  }
}
