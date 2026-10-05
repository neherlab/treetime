#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::project_root;
  use approx::assert_relative_eq;
  use pretty_assertions::assert_eq;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::io::json::json_read_str;
  use util_augur_node_data_json::AugurNodeDataJsonRefine;

  #[test]
  fn test_augur_node_data_optimize_branch_lengths() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");

    assert_relative_eq!(data.nodes["leaf_a"].branch_length, 0.005, max_relative = 1e-10);
    assert_relative_eq!(data.nodes["leaf_b"].branch_length, 0.010, max_relative = 1e-10);
    assert_relative_eq!(data.nodes["root"].branch_length, 0.0);
  }

  #[test]
  fn test_augur_node_data_optimize_only_branch_length_field() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");

    for (name, node) in &data.nodes {
      assert!(node.confidence.is_none(), "{name}: confidence must be omitted");
      assert!(node.numdate.is_none(), "{name}: numdate must be omitted");
      assert!(node.clock_length.is_none(), "{name}: clock_length must be omitted");
      assert!(
        node.mutation_length.is_none(),
        "{name}: mutation_length must be omitted"
      );
      assert!(node.raw_date.is_none(), "{name}: raw_date must be omitted");
      assert!(node.date.is_none(), "{name}: date must be omitted");
      assert!(node.date_inferred.is_none(), "{name}: date_inferred must be omitted");
      assert!(
        node.num_date_confidence.is_none(),
        "{name}: num_date_confidence must be omitted"
      );
    }
  }

  #[test]
  fn test_augur_node_data_optimize_metadata() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");

    assert!(data.metadata.clock.is_none());
    assert_eq!(data.metadata.alignment.as_deref(), Some("aln.fasta"));
    assert_eq!(data.metadata.input_tree.as_deref(), Some("tree.nwk"));
  }

  #[test]
  fn test_augur_node_data_optimize_generated_by() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");
    let generated_by = data.generated_by.unwrap();
    assert_eq!(generated_by.program, "treetime");
    assert_eq!(generated_by.version, env!("CARGO_PKG_VERSION"));
  }

  #[test]
  fn test_augur_node_data_optimize_roundtrip() {
    let json_str = helpers::write_json("(leaf_a:0.005,leaf_b:0.010)root;");

    let original: serde_json::Value = serde_json::from_str(&json_str).unwrap();
    let typed: AugurNodeDataJsonRefine = json_read_str(&json_str).unwrap();
    let roundtripped: serde_json::Value = serde_json::to_value(&typed).unwrap();

    assert_eq!(original, roundtripped);
  }

  #[test]
  fn test_augur_node_data_optimize_confidence_from_float_label() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)0.999:0.003;");

    assert_eq!(
      data.nodes["NODE_0000000"].confidence,
      Some(0.999),
      "Internal node with float label should emit confidence"
    );
    assert!(
      data.nodes["leaf_a"].confidence.is_none(),
      "Leaf nodes should not have confidence from float labels"
    );
    assert!(
      data.nodes["leaf_b"].confidence.is_none(),
      "Leaf nodes should not have confidence from float labels"
    );
  }

  #[test]
  fn test_augur_node_data_optimize_no_confidence_for_text_label() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");

    assert!(
      data.nodes["root"].confidence.is_none(),
      "Text-labeled internal node should not have confidence"
    );
  }

  #[test]
  fn test_augur_node_data_optimize_mutations_mode_branch_length_is_count() {
    let data = helpers::write_and_read_with_mutations("(leaf_a:0.005,leaf_b:0.010)root;", &[(0, 3), (1, 7)]);

    assert_relative_eq!(data.nodes["leaf_a"].branch_length, 3.0);
    assert_relative_eq!(data.nodes["leaf_b"].branch_length, 7.0);
    assert_relative_eq!(data.nodes["root"].branch_length, 0.0);
  }

  #[test]
  fn test_augur_node_data_optimize_mutations_mode_no_mutation_length() {
    let data = helpers::write_and_read_with_mutations("(leaf_a:0.005,leaf_b:0.010)root;", &[(0, 3), (1, 7)]);

    for node in data.nodes.values() {
      assert!(node.mutation_length.is_none());
    }
  }

  #[test]
  fn test_augur_node_data_optimize_end_to_end() {
    use treetime::alphabet::alphabet::Alphabet;
    use treetime::cancel::NoopCancel;
    use treetime::gtr::get_gtr::GtrModelName;
    use treetime::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
    use treetime::optimize::pipeline::{self, OptimizeInput, OptimizeParams};
    use treetime::progress::NoopProgress;
    use treetime_io::fasta::fasta_read_file;
    use treetime_io::nwk::nwk_read_file;

    let root = project_root();
    let alphabet = Alphabet::default();
    let nwk_parsed = nwk_read_file(root.join("data/flu/h3n2/20/tree.nwk")).unwrap();
    let confidences = nwk_parsed.confidences();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let sequences: Vec<AlignmentRecord> = fasta_read_file(root.join("data/flu/h3n2/20/aln.fasta.xz"), &alphabet)
      .unwrap()
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let params = OptimizeParams {
      model: GtrModelName::default(),
      dense: None,
      max_iter: 2,
      dp: 0.1,
      damping: 0.75,
      opt_method: BranchOptMethod::default(),
      initial_guess: InitialGuessMode::default(),
      no_indels: false,
      reroot_spec: None,
      topology_ops: TopologyOps::default(),
    };
    let input = OptimizeInput {
      graph,
      alphabet,
      sequences,
      branch_lengths,
    };

    let output = pipeline::run(&params, input, &names, &NoopCancel, &NoopProgress, &NoopProgress).unwrap();

    let data = helpers::build_augur_node_data_json_from_output(
      &output,
      &confidences,
      Some(std::path::Path::new("aln.fasta")),
      Some(std::path::Path::new("tree.nwk")),
    );

    assert!(!data.nodes.is_empty(), "node data must contain nodes");
    for (name, node) in &data.nodes {
      assert!(
        node.branch_length.is_finite() && node.branch_length >= 0.0,
        "{name}: branch_length must be finite and non-negative, got {}",
        node.branch_length
      );
    }
  }

  mod helpers {
    use app_output::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, Divergence, TreeSequences};
    use app_output::augur_node_data_refine::{RefineRun, build_augur_node_data_refine};
    use std::collections::BTreeMap;
    use std::path::Path;
    use treetime::optimize::pipeline::OptimizeOutput;
    use treetime::seq::mutation::Mutation;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::Seq;
    use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
    use util_augur_node_data_json::AugurNodeDataJsonRefine;

    pub(super) fn write_json(nwk: &str) -> String {
      let parse = nwk_read(nwk.as_bytes()).unwrap();
      let data = build(
        &parse.graph,
        &parse.names(),
        &parse.branch_lengths,
        &parse.confidences(),
        None,
        Some(Path::new("aln.fasta")),
        Some(Path::new("tree.nwk")),
      );
      json_write_str(&data, JsonPretty(true)).unwrap()
    }

    pub(super) fn write_and_read(nwk: &str) -> AugurNodeDataJsonRefine {
      json_read_str(write_json(nwk)).unwrap()
    }

    pub(super) fn write_and_read_with_mutations(nwk: &str, edge_counts: &[(usize, usize)]) -> AugurNodeDataJsonRefine {
      let parse = nwk_read(nwk.as_bytes()).unwrap();
      let edges = parse.graph.get_edges().collect::<Vec<_>>();
      let counts: BTreeMap<GraphEdgeKey, usize> = edge_counts
        .iter()
        .map(|&(idx, count)| (edges[idx].key(), count))
        .collect();
      let data = build(
        &parse.graph,
        &parse.names(),
        &parse.branch_lengths,
        &parse.confidences(),
        Some(&counts),
        Some(Path::new("aln.fasta")),
        Some(Path::new("tree.nwk")),
      );
      json_read_str(json_write_str(&data, JsonPretty(true)).unwrap()).unwrap()
    }

    pub(super) fn build_augur_node_data_json_from_output(
      output: &OptimizeOutput,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
      alignment: Option<&Path>,
      input_tree: Option<&Path>,
    ) -> AugurNodeDataJsonRefine {
      let data = build(
        &output.graph,
        &output.names,
        &output.branch_lengths,
        confidences,
        None,
        alignment,
        input_tree,
      );
      json_read_str(json_write_str(&data, JsonPretty(true)).unwrap()).unwrap()
    }

    fn build(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
      confidences: &BTreeMap<GraphNodeKey, Option<f64>>,
      mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
      alignment: Option<&Path>,
      input_tree: Option<&Path>,
    ) -> AugurNodeDataJsonRefine {
      let root_sequence = Seq::try_from_str("A").unwrap();
      let edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>> =
        graph.get_edges().map(|edge| (edge.key(), vec![])).collect();
      let annotated = AnnotatedGraph {
        graph,
        names,
        divergence_branch_lengths: branch_lengths,
        time_branch_lengths: None,
        divergence: Divergence::CumulativeBranchLength,
        branch_support: Some(confidences),
        sequences: mutation_counts.map(|mutation_counts| TreeSequences {
          root_sequence: &root_sequence,
          edge_mutations: &edge_mutations,
          mutation_counts: Some(mutation_counts),
          amino_acids: None,
        }),
        dates: None,
        traits: None,
      };
      let run = RefineRun {
        alignment,
        input_tree,
        clock_model: None,
        branch_support: Some(confidences),
      };
      build_augur_node_data_refine(&AnnotatedTreeView::new(&annotated).unwrap(), &run)
    }
  }
}
