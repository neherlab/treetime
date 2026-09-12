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
    assert_eq!(checked, 3); // leaf_a, leaf_b, and root (root's clock_length is 0.0)
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
  fn test_augur_node_data_timetree_root_has_zero_branch_fields() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

    // Root has no parent edge, so branch_length, clock_length, and mutation_length
    // are all zero. augur `export v2` requires mutation_length on every node
    // (including the root) to compute divergence. Internal/root nodes are always
    // date_inferred.
    let root = &data.nodes["root"];
    assert_relative_eq!(root.branch_length, 0.0);
    assert_relative_eq!(root.clock_length.unwrap(), 0.0);
    assert_relative_eq!(root.mutation_length.unwrap(), 0.0);
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
  fn test_augur_node_data_timetree_mutations_mode_root_has_zero_mutation_length() {
    let case = helpers::sample_case();
    let data = case.write_and_read_with_mutations(&[(0, 3), (1, 7)]);

    // The root has no parent edge, so its mutation_length is zero regardless of
    // the mutations-mode counts, which apply only to child edges.
    assert_relative_eq!(data.nodes["root"].mutation_length.unwrap(), 0.0);
    assert_relative_eq!(data.nodes["root"].branch_length, 0.0);
  }

  mod helpers {
    use crate::clock::clock_model::{ClockModel, ClockModelStats, RegressionStats};
    use crate::commands::timetree::output::augur_node_data::build_augur_node_data_json;
    use crate::commands::timetree::result::{TimetreeEdgeOut, TimetreeNodeOut};
    use crate::timetree::confidence::NodeConfidenceInterval;
    use ndarray::array;
    use std::collections::BTreeMap;
    use std::path::Path;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::dates_csv::{DateConstraint, DateRange, DateValue, DatesMap};
    use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
    use util_augur_node_data_json::AugurNodeDataJsonRefine;

    pub struct SampleCase {
      pub graph: Graph,
      pub names: BTreeMap<GraphNodeKey, Option<String>>,
      pub times: BTreeMap<GraphNodeKey, Option<f64>>,
      pub clock_model: ClockModel,
      pub dates: DatesMap,
      pub intervals: Vec<NodeConfidenceInterval>,
      pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
    }

    impl SampleCase {
      pub fn write_json(&self) -> String {
        let data = build_augur_node_data_json(
          &self.graph,
          &timetree_nodes(&self.graph, &self.names, &self.times),
          &timetree_edges(&self.graph, &self.branch_lengths),
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
          &timetree_nodes(&self.graph, &self.names, &self.times),
          &timetree_edges(&self.graph, &self.branch_lengths),
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
      let mut graph = Graph::new();
      let mut names = BTreeMap::new();
      let mut times = BTreeMap::new();
      let root_key = graph.add_node();
      names.insert(root_key, Some("root".to_owned()));
      times.insert(root_key, Some(2000.0));
      let leaf_a_key = graph.add_node();
      names.insert(leaf_a_key, Some("leaf_a".to_owned()));
      times.insert(leaf_a_key, Some(2005.0));
      let leaf_b_key = graph.add_node();
      names.insert(leaf_b_key, Some("leaf_b".to_owned()));
      times.insert(leaf_b_key, Some(2010.0));
      let edge_a_key = graph.add_edge(root_key, leaf_a_key).unwrap();
      let edge_b_key = graph.add_edge(root_key, leaf_b_key).unwrap();
      graph.build().unwrap();

      let branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = maplit::btreemap! {
        edge_a_key => Some(0.005),
        edge_b_key => Some(0.010),
      };

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
        names,
        times,
        clock_model,
        dates,
        intervals,
        branch_lengths,
      }
    }

    fn timetree_nodes<D: Send + Sync>(
      graph: &Graph<D>,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      times: &BTreeMap<GraphNodeKey, Option<f64>>,
    ) -> BTreeMap<GraphNodeKey, TimetreeNodeOut> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let key = node.read_arc().key();
          (
            key,
            TimetreeNodeOut {
              name: names.get(&key).cloned().flatten(),
              desc: None,
              // This fixture builds its graph node by node with no input-tree branch support, so
              // production's parse-time confidence map would surface None for every node here too.
              confidence: None,
              time: times.get(&key).copied().flatten(),
              div: 0.0,
              is_outlier: false,
              bad_branch: false,
              // Rate-susceptibility dates are produced only by the confidence pass and threaded as a
              // value map; this fixture graph runs no such pass, so production surfaces None here too.
              rate_susceptibility_dates: None,
            },
          )
        })
        .collect()
    }

    fn timetree_edges<D: Send + Sync>(
      graph: &Graph<D>,
      branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    ) -> BTreeMap<GraphEdgeKey, TimetreeEdgeOut> {
      graph
        .get_edges()
        .iter()
        .map(|edge| {
          let key = edge.read_arc().key();
          (
            key,
            TimetreeEdgeOut {
              branch_length: branch_lengths[&key],
              time_length: None,
              clock_branch_length: None,
              // Strict-clock test graph: the relaxed-clock multiplier is its default 1.0, matching
              // what production reads from the threaded edge state.
              gamma: 1.0,
            },
          )
        })
        .collect()
    }
  }
}
