#[cfg(test)]
mod tests {
  use approx::assert_relative_eq;
  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::json_read_str;
  use util_augur_node_data_json::AugurNodeDataJsonRefine;

  // --- Full output: per-node dates, branch metrics, confidence, and clock block ---
  //
  // Field semantics match augur refine (traced against TreeTime's final state):
  //   mutation_length = ML divergence length (subs/site) = parent edge base length
  //   branch_length   = clock_length = child.numdate - parent.numdate (years)

  #[test]
  fn test_augur_node_data_timetree_full_output() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

    // leaf_a: edge divergence 0.005, 5-year branch, exact date -> not inferred.
    let leaf_a = &data.nodes["leaf_a"];
    assert_relative_eq!(leaf_a.mutation_length.unwrap(), 0.005); // ML divergence (edge base length)
    assert_relative_eq!(leaf_a.branch_length, 5.0); // child.time - parent.time = 2005 - 2000
    assert_relative_eq!(leaf_a.clock_length.unwrap(), 5.0); // identical to branch_length
    assert_relative_eq!(leaf_a.numdate.unwrap(), 2005.0);
    assert_eq!(leaf_a.date.as_deref(), Some("2005-01-01"));
    assert_eq!(leaf_a.raw_date.as_deref(), Some("2005"));
    assert_eq!(leaf_a.date_inferred, Some(false));
    let ci = leaf_a.num_date_confidence.unwrap();
    assert_relative_eq!(ci[0], 2004.0);
    assert_relative_eq!(ci[1], 2006.0);

    // leaf_b: edge divergence 0.010, 10-year branch, uncertain date -> inferred.
    let leaf_b = &data.nodes["leaf_b"];
    assert_relative_eq!(leaf_b.mutation_length.unwrap(), 0.010);
    assert_relative_eq!(leaf_b.branch_length, 10.0);
    assert_relative_eq!(leaf_b.clock_length.unwrap(), 10.0);
    assert_eq!(leaf_b.raw_date.as_deref(), Some("2010-XX-XX"));
    assert_eq!(leaf_b.date_inferred, Some(true));
    assert!(leaf_b.num_date_confidence.is_none());

    // clock block from the regression covariance.
    let clock = data.metadata.clock.as_ref().unwrap();
    assert_relative_eq!(clock.rate, 0.002);
    assert_relative_eq!(clock.intercept, -4.0);
    assert_relative_eq!(clock.rtt_tmrca, 2000.0); // -(-4.0) / 0.002
    let cov = clock.cov.as_ref().unwrap();
    assert_eq!(cov, &vec![vec![1e-8, 0.0], vec![0.0, 0.5]]);
    assert_relative_eq!(clock.rate_std.unwrap(), 1e-4); // sqrt(cov[0,0])

    assert_eq!(data.metadata.input_tree.as_deref(), Some("tree.nwk"));
    assert_eq!(data.metadata.alignment.as_deref(), Some("aln.fasta"));
  }

  // --- Mapping invariants that must never silently regress (augur parity) ---

  #[test]
  fn test_augur_node_data_timetree_branch_length_equals_clock_length() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

    // augur sets node.branch_length = node.clock_length (clock_tree.py:925). For
    // every node with a clock_length the two fields are identical. Both leaves
    // carry distinct non-zero values, so this is not a vacuous 0.0 == 0.0 check.
    let mut checked = 0;
    for node in data.nodes.values() {
      if let Some(clock_length) = node.clock_length {
        assert_relative_eq!(node.branch_length, clock_length);
        checked += 1;
      }
    }
    assert_eq!(checked, 2); // leaf_a and leaf_b (root has no clock_length)
  }

  #[test]
  fn test_augur_node_data_timetree_mutation_length_is_divergence() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

    // mutation_length is the ML divergence (parent edge base length), distinct
    // from the time-valued branch_length. Asserting inequality guards against a
    // regression to emitting the time value (or count/site) here.
    let leaf_a = &data.nodes["leaf_a"];
    assert_relative_eq!(leaf_a.mutation_length.unwrap(), 0.005);
    assert!((leaf_a.mutation_length.unwrap() - leaf_a.branch_length).abs() > 1.0);
  }

  #[test]
  fn test_augur_node_data_timetree_root_has_no_branch_fields() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

    // Root has no parent edge: branch_length defaults to 0.0, clock_length and
    // mutation_length are omitted. Internal/root nodes are always date_inferred.
    let root = &data.nodes["root"];
    assert_relative_eq!(root.branch_length, 0.0);
    assert!(root.clock_length.is_none());
    assert!(root.mutation_length.is_none());
    assert!(root.raw_date.is_none());
    assert_eq!(root.date_inferred, Some(true));
    assert_relative_eq!(root.numdate.unwrap(), 2000.0);
  }

  #[test]
  fn test_augur_node_data_timetree_roundtrip() {
    let case = helpers::sample_case();
    let json_str = case.write_json();

    let original: serde_json::Value = serde_json::from_str(&json_str).unwrap();
    let typed: AugurNodeDataJsonRefine = json_read_str(&json_str).unwrap();
    let roundtripped: serde_json::Value = serde_json::to_value(&typed).unwrap();

    assert_eq!(original, roundtripped);
  }

  #[test]
  fn test_augur_node_data_timetree_generated_by() {
    let case = helpers::sample_case();
    let data = case.write_and_read();
    let generated_by = data.generated_by.unwrap();
    assert_eq!(generated_by.program, "treetime");
    assert_eq!(generated_by.version, env!("CARGO_PKG_VERSION"));
  }

  // --- Divergence units: mutations mode ---

  #[test]
  fn test_augur_node_data_timetree_mutations_mode_mutation_length_is_count() {
    let case = helpers::sample_case();
    let data = case.write_and_read_with_mutations(&[(0, 3), (1, 7)]);

    assert_relative_eq!(data.nodes["leaf_a"].mutation_length.unwrap(), 3.0);
    assert_relative_eq!(data.nodes["leaf_b"].mutation_length.unwrap(), 7.0);
  }

  #[test]
  fn test_augur_node_data_timetree_mutations_mode_branch_length_stays_years() {
    let case = helpers::sample_case();
    let data = case.write_and_read_with_mutations(&[(0, 3), (1, 7)]);

    assert_relative_eq!(data.nodes["leaf_a"].branch_length, 5.0);
    assert_relative_eq!(data.nodes["leaf_b"].branch_length, 10.0);
  }

  #[test]
  fn test_augur_node_data_timetree_mutations_mode_root_has_no_mutation_length() {
    let case = helpers::sample_case();
    let data = case.write_and_read_with_mutations(&[(0, 3), (1, 7)]);

    assert!(data.nodes["root"].mutation_length.is_none());
    assert_relative_eq!(data.nodes["root"].branch_length, 0.0);
  }

  mod helpers {
    use crate::clock::clock_model::{ClockModel, ClockModelStats, RegressionStats};
    use crate::commands::timetree::output::augur_node_data::build_augur_node_data_json;
    use crate::partition::timetree::GraphTimetree;
    use crate::payload::timetree::{EdgeTimetree, NodeTimetree};
    use crate::timetree::confidence::NodeConfidenceInterval;
    use ndarray::array;
    use std::collections::BTreeMap;
    use std::path::Path;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::node::Named;
    use treetime_io::dates_csv::{DateConstraint, DateRange, DateValue, DatesMap};
    use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
    use util_augur_node_data_json::AugurNodeDataJsonRefine;

    pub struct SampleCase {
      pub graph: GraphTimetree,
      pub clock_model: ClockModel,
      pub dates: DatesMap,
      pub intervals: Vec<NodeConfidenceInterval>,
    }

    impl SampleCase {
      pub fn write_json(&self) -> String {
        let data = build_augur_node_data_json(
          &self.graph,
          &self.clock_model,
          Some(&self.intervals),
          Some(&self.dates),
          Some(Path::new("aln.fasta")),
          Some(Path::new("tree.nwk")),
          None,
        )
        .unwrap();
        json_write_str(&data, JsonPretty(true)).unwrap()
      }

      pub fn write_and_read(&self) -> AugurNodeDataJsonRefine {
        json_read_str(self.write_json()).unwrap()
      }

      pub fn write_and_read_with_mutations(&self, edge_counts: &[(usize, usize)]) -> AugurNodeDataJsonRefine {
        let edges = self.graph.get_edges();
        let counts: BTreeMap<GraphEdgeKey, usize> = edge_counts
          .iter()
          .map(|&(idx, count)| (edges[idx].read_arc().key(), count))
          .collect();
        let data = build_augur_node_data_json(
          &self.graph,
          &self.clock_model,
          Some(&self.intervals),
          Some(&self.dates),
          Some(Path::new("aln.fasta")),
          Some(Path::new("tree.nwk")),
          Some(&counts),
        )
        .unwrap();
        json_read_str(json_write_str(&data, JsonPretty(true)).unwrap()).unwrap()
      }
    }

    /// Build a 3-node tree (root -> leaf_a, root -> leaf_b) with node times, edge
    /// divergence lengths (the `mutation_length` source), a regression clock model,
    /// date constraints (exact for leaf_a, uncertain for leaf_b), and CIs for root
    /// and leaf_a.
    pub fn sample_case() -> SampleCase {
      let mut graph = GraphTimetree::new();
      let root_key = graph.add_node(make_node("root", 2000.0));
      let leaf_a_key = graph.add_node(make_node("leaf_a", 2005.0));
      let leaf_b_key = graph.add_node(make_node("leaf_b", 2010.0));
      graph.add_edge(root_key, leaf_a_key, make_edge(0.005)).unwrap();
      graph.add_edge(root_key, leaf_b_key, make_edge(0.010)).unwrap();
      graph.build().unwrap();

      // Regression covariance over [rate, intercept]: cov[0,0] is the rate variance.
      let stats = ClockModelStats::Estimated(RegressionStats {
        chisq: 0.0,
        r_val: 0.99,
        hessian: array![[0.0, 0.0], [0.0, 0.0]],
        cov: array![[1e-8, 0.0], [0.0, 0.5]],
      });
      let clock_model = ClockModel::for_testing_with_stats(0.002, -4.0, stats);

      let dates: DatesMap = maplit::btreemap! {
        "leaf_a".to_owned() => Some(DateConstraint::exact(2005.0)),
        "leaf_b".to_owned() => Some(DateConstraint {
          raw: "2010-XX-XX".to_owned(),
          value: DateValue::Uncertain(DateRange { start: 2010.0, end: 2011.0 }),
        }),
      };

      let intervals = vec![
        NodeConfidenceInterval {
          key: root_key,
          name: "root".to_owned(),
          date: 2000.0,
          lower: 1998.0,
          upper: 2002.0,
        },
        NodeConfidenceInterval {
          key: leaf_a_key,
          name: "leaf_a".to_owned(),
          date: 2005.0,
          lower: 2004.0,
          upper: 2006.0,
        },
      ];

      SampleCase {
        graph,
        clock_model,
        dates,
        intervals,
      }
    }

    fn make_node(name: &str, time: f64) -> NodeTimetree {
      let mut node = NodeTimetree::default();
      node.base.set_name(Some(name));
      node.time = Some(time);
      node
    }

    fn make_edge(branch_length: f64) -> EdgeTimetree {
      let mut edge = EdgeTimetree::default();
      edge.base.branch_length = Some(branch_length);
      edge
    }
  }
}
