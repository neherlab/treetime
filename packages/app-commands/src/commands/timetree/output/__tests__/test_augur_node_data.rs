#[cfg(test)]
mod tests {
  use approx::assert_relative_eq;
  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::json_read_str;
  use util_augur_node_data_json::AugurNodeDataJsonRefine;

  #[test]
  fn test_augur_node_data_timetree_full_output() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

    let leaf_a = &data.nodes["leaf_a"];
    assert_relative_eq!(leaf_a.mutation_length.unwrap(), 0.005);
    assert_relative_eq!(leaf_a.branch_length, 5.0);
    assert_relative_eq!(leaf_a.clock_length.unwrap(), 5.0);
    assert_relative_eq!(leaf_a.numdate.unwrap(), 2005.0);
    assert_eq!(leaf_a.date.as_deref(), Some("2005-01-01"));
    assert_eq!(leaf_a.raw_date.as_deref(), Some("2005"));
    assert_eq!(leaf_a.date_inferred, Some(false));
    let ci = leaf_a.num_date_confidence.unwrap();
    assert_relative_eq!(ci[0], 2004.0);
    assert_relative_eq!(ci[1], 2006.0);

    let leaf_b = &data.nodes["leaf_b"];
    assert_relative_eq!(leaf_b.mutation_length.unwrap(), 0.010);
    assert_relative_eq!(leaf_b.branch_length, 10.0);
    assert_relative_eq!(leaf_b.clock_length.unwrap(), 10.0);
    assert_eq!(leaf_b.raw_date.as_deref(), Some("2010-XX-XX"));
    assert_eq!(leaf_b.date_inferred, Some(true));
    assert!(leaf_b.num_date_confidence.is_none());

    let clock = data.metadata.clock.as_ref().unwrap();
    assert_relative_eq!(clock.rate, 0.002);
    assert_relative_eq!(clock.intercept, -4.0);
    assert_relative_eq!(clock.rtt_tmrca, 2000.0);
    let cov = clock.cov.as_ref().unwrap();
    assert_eq!(cov, &vec![vec![1e-8, 0.0], vec![0.0, 0.5]]);
    assert_relative_eq!(clock.rate_std.unwrap(), 1e-4);

    assert_eq!(data.metadata.input_tree.as_deref(), Some("tree.nwk"));
    assert_eq!(data.metadata.alignment.as_deref(), Some("aln.fasta"));
  }

  #[test]
  fn test_augur_node_data_timetree_branch_length_equals_clock_length() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

    let mut checked = 0;
    for node in data.nodes.values() {
      if let Some(clock_length) = node.clock_length {
        assert_relative_eq!(node.branch_length, clock_length);
        checked += 1;
      }
    }
    assert_eq!(checked, 3);
  }

  #[test]
  fn test_augur_node_data_timetree_mutation_length_is_divergence() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

    let leaf_a = &data.nodes["leaf_a"];
    assert_relative_eq!(leaf_a.mutation_length.unwrap(), 0.005);
    assert!((leaf_a.mutation_length.unwrap() - leaf_a.branch_length).abs() > 1.0);
  }

  #[test]
  fn test_augur_node_data_timetree_root_has_zero_branch_fields() {
    let case = helpers::sample_case();
    let data = case.write_and_read();

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

    assert_relative_eq!(data.nodes["root"].mutation_length.unwrap(), 0.0);
    assert_relative_eq!(data.nodes["root"].branch_length, 0.0);
  }

  mod helpers {
    use crate::__tests__::test_support::tests::NUC_ALPHABET;
    use app_output::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, Divergence, TreeDates, TreeSequences};
    use app_output::augur_node_data_refine::{RefineRun, build_augur_node_data_refine};
    use indoc::indoc;
    use std::collections::{BTreeMap, BTreeSet};
    use std::path::Path;
    use treetime::clock::clock_model::{ClockModel, ClockRegression};
    use treetime::seq::mutation::Mutation;
    use treetime::timetree::confidence::NodeConfidenceInterval;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::dates_csv::{DateConstraint, DateRange, DateValue};
    use treetime_primitives::Seq;
    use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
    use util_augur_node_data_json::AugurNodeDataJsonRefine;

    pub(super) struct SampleCase {
      graph: Graph,
      names: BTreeMap<GraphNodeKey, Option<String>>,
      times: BTreeMap<GraphNodeKey, Option<f64>>,
      clock_model: ClockModel,
      dates: BTreeMap<GraphNodeKey, DateConstraint>,
      intervals: Vec<NodeConfidenceInterval>,
      branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
    }

    impl SampleCase {
      pub(crate) fn write_json(&self) -> String {
        json_write_str(&self.build(None), JsonPretty(true)).unwrap()
      }

      pub(crate) fn write_and_read(&self) -> AugurNodeDataJsonRefine {
        json_read_str(self.write_json()).unwrap()
      }

      pub(crate) fn write_and_read_with_mutations(&self, edge_counts: &[(usize, usize)]) -> AugurNodeDataJsonRefine {
        let edges = self.graph.get_edges().collect::<Vec<_>>();
        let counts: BTreeMap<GraphEdgeKey, usize> = edge_counts
          .iter()
          .map(|&(idx, count)| (edges[idx].key(), count))
          .collect();
        json_read_str(json_write_str(&self.build(Some(&counts)), JsonPretty(true)).unwrap()).unwrap()
      }

      fn build(&self, mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>) -> AugurNodeDataJsonRefine {
        let time_lengths = time_lengths(&self.graph, &self.times);
        let confidence: BTreeMap<GraphNodeKey, [f64; 2]> = self
          .intervals
          .iter()
          .map(|interval| (interval.key, [interval.lower, interval.upper]))
          .collect();
        let excluded = BTreeSet::new();
        let root_sequence = Seq::try_from_str("A").unwrap();
        let edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>> =
          self.graph.get_edges().map(|edge| (edge.key(), vec![])).collect();
        let annotated = AnnotatedGraph {
          graph: &self.graph,
          names: &self.names,
          divergence_branch_lengths: &self.branch_lengths,
          time_branch_lengths: Some(&time_lengths),
          divergence: Divergence::CumulativeBranchLength,
          sequences: mutation_counts.map(|mutation_counts| TreeSequences {
            alphabet: &NUC_ALPHABET,
            root_sequence: &root_sequence,
            edge_mutations: &edge_mutations,
            mutation_counts: Some(mutation_counts),
            amino_acids: None,
          }),
          dates: Some(TreeDates {
            num_date: &self.times,
            confidence: Some(&confidence),
            excluded: &excluded,
            input_dates: Some(&self.dates),
          }),
          traits: None,
        };
        let run = RefineRun {
          alignment: Some(Path::new("aln.fasta")),
          input_tree: Some(Path::new("tree.nwk")),
          clock_model: Some(&self.clock_model),
        };
        build_augur_node_data_refine(&AnnotatedTreeView::new(&annotated).unwrap(), &run)
      }
    }

    pub(super) fn sample_case() -> SampleCase {
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

      let clock_model = sample_clock_model();

      let dates: BTreeMap<GraphNodeKey, DateConstraint> = maplit::btreemap! {
        leaf_a_key => DateConstraint::exact(2005.0),
        leaf_b_key => DateConstraint {
          raw: "2010-XX-XX".to_owned(),
          value: DateValue::Uncertain(DateRange { start: 2010.0, end: 2011.0 }),
        },
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

    fn sample_clock_model() -> ClockModel {
      let regression: ClockRegression = json_read_str(indoc! {r#"{
        "clock_rate": 0.002,
        "intercept": -4.0,
        "chisq": 0.0,
        "r_val": 0.99,
        "hessian": [[0.0, 0.0], [0.0, 0.0]],
        "cov": [[1e-8, 0.0], [0.0, 0.5]]
      }"#})
      .unwrap();
      ClockModel::from_regression(&regression).unwrap()
    }

    fn time_lengths(graph: &Graph, times: &BTreeMap<GraphNodeKey, Option<f64>>) -> BTreeMap<GraphEdgeKey, Option<f64>> {
      graph
        .get_edges()
        .map(|edge| {
          let parent = times.get(&edge.source()).copied().flatten();
          let child = times.get(&edge.target()).copied().flatten();
          (edge.key(), parent.zip(child).map(|(parent, child)| child - parent))
        })
        .collect()
    }
  }
}
