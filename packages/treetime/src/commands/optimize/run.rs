use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::compress_sequences;
use crate::commands::ancestral::marginal_dense::run_marginal_dense;
use crate::commands::ancestral::marginal_sparse::run_marginal_sparse;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::optimize_dense::run_optimize_dense;
use crate::commands::optimize::optimize_sparse::run_optimize_sparse;
use crate::graph::edge::{GraphEdge, NumMuts, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use crate::gtr::get_gtr::{get_gtr, get_gtr_dense};
use crate::io::fasta::read_many_fasta;
use crate::io::nex::{NexWriteOptions, nex_write_file};
use crate::io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use crate::make_error;
use crate::representation::graph_dense::DenseGraph;
use crate::representation::graph_sparse::SparseGraph;
use crate::representation::infer_dense::infer_dense;
use crate::representation::partitions_likelihood::{PartitionLikelihood, PartitionLikelihoodWithAln};
use crate::representation::partitions_parsimony::PartitionParsimonyWithAln;
use crate::utils::float_fmt::float_to_significant_digits;
use eyre::Report;
use itertools::Itertools;
use log::debug;
use serde::Serialize;
use std::path::Path;

// The initial guess for dense is not working well, but optimization works without revisit after settling on optimization algorithm
// use super::optimize_dense::initial_guess;
use super::optimize_sparse::initial_guess_sparse;

#[derive(Clone, Debug, Default)]
pub struct TreetimeOptimizeParams {
  pub sample_from_profile: bool,
  pub fixed_pi: bool,
}

fn validate_args(args: &TreetimeOptimizeArgs) -> Result<(), Report> {
  if args.prune_empty && args.input_fastas.is_empty() {
    return make_error!(
      "The --prune-empty requires --aln. Without sequence data, it's not possible to determine which branches lack mutations."
    );
  }

  Ok(())
}

pub fn run_optimize(args: &TreetimeOptimizeArgs) -> Result<(), Report> {
  validate_args(args)?;

  let TreetimeOptimizeArgs {
    input_fastas,
    tree,
    alphabet,
    model_name,
    dense,
    outdir,
    max_iter,
    dp,
    prune_short,
    prune_empty,
  } = args;

  let dense = dense.unwrap_or_else(infer_dense);

  let treat_gap_as_unknown = dense;
  let alphabet = Alphabet::new(alphabet.unwrap_or_default(), treat_gap_as_unknown)?;

  // TODO: avoid reading all sequences into memory somehow?
  let aln = read_many_fasta(input_fastas, &alphabet)?;

  // TODO: refactor to reduce duplication with `ancestral` as well as within the branches of this conditional
  if !dense {
    let mut graph: SparseGraph = nwk_read_file(tree)?;
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet.clone(), aln)?];
    let partitions = compress_sequences(&graph, partitions)?;

    let gtr = get_gtr(model_name, &alphabet, &graph)?;
    let partitions = partitions
          .into_iter()
          .map(|part| PartitionLikelihood::from_parsimony(gtr.clone(), part)) // FIXME: avoid cloning
          .collect_vec();

    initial_guess_sparse(&graph, &partitions);
    let mut lh_prev = f64::MIN;
    for i in 0..*max_iter {
      let lh = run_marginal_sparse(&graph, &partitions)?;
      debug!("Iteration {}: likelihood {}", i + 1, float_to_significant_digits(lh, 7));
      if (lh - lh_prev).abs() < dp.abs() {
        break;
      }
      run_optimize_sparse(&graph, &partitions)?;
      lh_prev = lh;
    }

    collapse_short_edges(&mut graph, *prune_short, *prune_empty)?;

    write_graph(outdir, &graph)?;
  } else {
    let mut graph: DenseGraph = nwk_read_file(tree)?;
    let gtr = get_gtr_dense(model_name, &alphabet, &graph)?;

    let partitions_waln = vec![PartitionLikelihoodWithAln::new(gtr, alphabet, aln)?];
    let partitions = partitions_waln
      .iter()
      .map(|part| PartitionLikelihood::from(part.clone()))
      .collect_vec();
    let mut lh_prev = f64::MIN;
    for i in 0..*max_iter {
      // FIXME: avoid assigning sequences to the graph in every iteration
      let lh = run_marginal_dense(&graph, partitions_waln.clone(), false)?; // FIXME: avoid cloning

      // somehow, the initial guess makes it worse...
      // if i == 0 {
      //   initial_guess(&graph, &partitions);
      // }

      debug!("Iteration {}: likelihood {}", i + 1, float_to_significant_digits(lh, 7));
      if (lh - lh_prev).abs() < dp.abs() {
        break;
      }
      run_optimize_dense(&graph, &partitions)?;
      lh_prev = lh;
    }

    collapse_short_edges(&mut graph, *prune_short, *prune_empty)?;

    write_graph(outdir, &graph)?;
  }

  Ok(())
}

/// Collapse edges with weights below the given threshold.
/// NOTE: Leaf nodes are excluded from collapsing.
pub fn collapse_short_edges<N, E, D>(
  graph: &mut Graph<N, E, D>,
  prune_short: Option<f64>,
  prune_empty: bool,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge + Weighted + NumMuts,
  D: Sync + Send,
{
  #[allow(clippy::needless_collect)]
  let edges_to_collapse: Vec<_> = graph
    .get_edges()
    .iter()
    .filter_map(|edge| {
      let edge = edge.read_arc();
      let weight = edge.payload().read_arc().weight();
      let target_is_leaf = graph.is_leaf(edge.target());
      let should_prune_short = matches!((prune_short, weight), (Some(threshold), Some(weight)) if weight < threshold);
      let should_prune_empty = prune_empty && edge.payload().read_arc().num_muts() == Some(0);
      ((should_prune_short || should_prune_empty) && !target_is_leaf).then(|| edge.key())
    })
    .collect();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing edge: {edge_key}");
    graph.collapse_edge(edge_key)
  })
}

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + Serialize,
  E: GraphEdge + EdgeToNwk + Serialize,
  D: Send + Sync + Default + Serialize,
{
  // json_write_file(
  //   outdir.as_ref().join("annotated_tree.graph.json"),
  //   &graph,
  //   JsonPretty(true),
  // )?;

  nwk_write_file(
    outdir.as_ref().join("annotated_tree.nwk"),
    graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write_file(
    outdir.as_ref().join("annotated_tree.nexus"),
    graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::edge::NumMuts;
  use crate::graph::graph::Graph;
  use crate::graph::graph::tests::TestNode;
  use crate::io::nwk::{nwk_read_str, nwk_write_str};
  use pretty_assertions::assert_eq;
  use serde::{Deserialize, Serialize};

  #[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
  struct TestEdgeWithMuts {
    weight: Option<f64>,
    num_muts: Option<usize>,
  }

  impl TestEdgeWithMuts {
    fn new(weight: Option<f64>, num_muts: Option<usize>) -> Self {
      Self { weight, num_muts }
    }
  }

  impl GraphEdge for TestEdgeWithMuts {}

  impl Weighted for TestEdgeWithMuts {
    fn weight(&self) -> Option<f64> {
      self.weight
    }

    fn set_weight(&mut self, weight: Option<f64>) {
      self.weight = weight;
    }
  }

  impl NumMuts for TestEdgeWithMuts {
    fn num_muts(&self) -> Option<usize> {
      self.num_muts
    }
  }

  impl crate::io::nwk::EdgeFromNwk for TestEdgeWithMuts {
    fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
      Ok(Self::new(weight, None))
    }
  }

  impl EdgeToNwk for TestEdgeWithMuts {
    fn nwk_weight(&self) -> Option<f64> {
      self.weight
    }
  }

  #[test]
  fn test_collapse_short_edges_basic() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> = nwk_read_str("(A:0.0,B:0.1)root;")?;
    collapse_short_edges(&mut graph, Some(0.0), false)?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_with_threshold() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> = nwk_read_str("(A:0.01,B:0.02,C:0.1)root;")?;
    collapse_short_edges(&mut graph, Some(0.05), false)?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.01,B:0.02,C:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_preserves_large_edges() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> = nwk_read_str("(A:0.1,B:0.2)root;")?;
    collapse_short_edges(&mut graph, Some(0.01), false)?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,B:0.2)root;");
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_empty_graph() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> = Graph::new();
    collapse_short_edges(&mut graph, Some(0.0), false)?;
    assert!(graph.get_nodes().is_empty());
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_handles_none_weights() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> = nwk_read_str("(A:0.0,B:0.1)root;")?;
    collapse_short_edges(&mut graph, Some(0.0), false)?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_preserves_terminal_nodes() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> = nwk_read_str("(A:0.00001,B:0.1)root;")?;
    collapse_short_edges(&mut graph, Some(0.001), false)?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.00001,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_complex_tree() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> =
      nwk_read_str("(((A:0,B:0.1)internal1:0.00002,(C:0.00003,D:0.1)internal2:0.1)internal3:0.00004,E:0.00005)root;")?;
    collapse_short_edges(&mut graph, Some(0.01), false)?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(E:0.00005,(C:0.00003,D:0.1)internal2:0.1,A:0,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_prune_empty_preserves_leaves() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdgeWithMuts, ()>::new();

    let root = graph.add_node(TestNode(Some("root".to_owned())));
    let a = graph.add_node(TestNode(Some("A".to_owned())));
    let b = graph.add_node(TestNode(Some("B".to_owned())));

    // Edge with no mutations to leaf A (should be preserved)
    let root_to_a = graph.add_edge(root, a, TestEdgeWithMuts::new(Some(0.1), Some(0)))?;
    // Edge with mutations to leaf B
    let root_to_b = graph.add_edge(root, b, TestEdgeWithMuts::new(Some(0.1), Some(2)))?;

    graph.build()?;

    collapse_short_edges(&mut graph, None, true)?;

    // Both leaves should be preserved even if edge to A has no mutations
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_prune_empty_internal_nodes() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdgeWithMuts, ()>::new();

    let root = graph.add_node(TestNode(Some("root".to_owned())));
    let internal = graph.add_node(TestNode(Some("internal".to_owned())));
    let a = graph.add_node(TestNode(Some("A".to_owned())));
    let b = graph.add_node(TestNode(Some("B".to_owned())));

    // Internal edge with no mutations
    let root_to_internal = graph.add_edge(root, internal, TestEdgeWithMuts::new(Some(0.1), Some(0)))?;
    // Leaf edges with mutations
    let internal_to_a = graph.add_edge(internal, a, TestEdgeWithMuts::new(Some(0.1), Some(1)))?;
    let internal_to_b = graph.add_edge(internal, b, TestEdgeWithMuts::new(Some(0.1), Some(2)))?;

    graph.build()?;

    collapse_short_edges(&mut graph, None, true)?;

    // Internal node should be collapsed, but leaves preserved
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_prune_empty_none_mutations() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdgeWithMuts, ()>::new();

    let root = graph.add_node(TestNode(Some("root".to_owned())));
    let internal = graph.add_node(TestNode(Some("internal".to_owned())));
    let a = graph.add_node(TestNode(Some("A".to_owned())));

    // Edge with None mutations (should not be pruned - None means unknown, not zero)
    let root_to_internal = graph.add_edge(root, internal, TestEdgeWithMuts::new(Some(0.1), None))?;
    let internal_to_a = graph.add_edge(internal, a, TestEdgeWithMuts::new(Some(0.1), Some(1)))?;

    graph.build()?;

    collapse_short_edges(&mut graph, None, true)?;

    // Internal node should be preserved when mutations is None (unknown)
    assert_eq!(graph.get_nodes().len(), 3); // root, internal, A
    assert_eq!(graph.get_edges().len(), 2); // root->internal, internal->A

    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_prune_empty_simple_leaf_case() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdgeWithMuts, ()>::new();

    let root = graph.add_node(TestNode(Some("root".to_owned())));
    let a = graph.add_node(TestNode(Some("A".to_owned())));

    // Single leaf edge with no mutations should be preserved
    let root_to_a = graph.add_edge(root, a, TestEdgeWithMuts::new(Some(0.1), Some(0)))?;

    graph.build()?;

    collapse_short_edges(&mut graph, None, true)?;

    // Leaf should be preserved even with no mutations
    assert_eq!(graph.get_nodes().len(), 2); // root and A
    assert_eq!(graph.get_edges().len(), 1); // root->A

    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_combined_prune_short_and_empty() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdgeWithMuts, ()>::new();

    let root = graph.add_node(TestNode(Some("root".to_owned())));
    let internal1 = graph.add_node(TestNode(Some("internal1".to_owned())));
    let internal2 = graph.add_node(TestNode(Some("internal2".to_owned())));
    let a = graph.add_node(TestNode(Some("A".to_owned())));
    let b = graph.add_node(TestNode(Some("B".to_owned())));

    // Short edge (should be pruned by threshold)
    let root_to_internal1 = graph.add_edge(root, internal1, TestEdgeWithMuts::new(Some(0.001), Some(1)))?;
    // Empty edge (should be pruned by empty check)
    let root_to_internal2 = graph.add_edge(root, internal2, TestEdgeWithMuts::new(Some(0.1), Some(0)))?;
    // Leaf edges (should be preserved)
    let internal1_to_a = graph.add_edge(internal1, a, TestEdgeWithMuts::new(Some(0.1), Some(1)))?;
    let internal2_to_b = graph.add_edge(internal2, b, TestEdgeWithMuts::new(Some(0.1), Some(2)))?;

    graph.build()?;

    collapse_short_edges(&mut graph, Some(0.01), true)?;

    // Both internal nodes should be collapsed, leaves preserved
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_prune_short_threshold_exact() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> = nwk_read_str("(A:0.05,B:0.05,C:0.051)root;")?;
    collapse_short_edges(&mut graph, Some(0.05), false)?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.05,B:0.05,C:0.051)root;");
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_prune_short_threshold_below() -> Result<(), Report> {
    let mut graph: Graph<TestNode, TestEdgeWithMuts, ()> = nwk_read_str("(A:0.049,B:0.05,C:0.051)root;")?;
    collapse_short_edges(&mut graph, Some(0.05), false)?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    // All edges remain because A, B, C are leaves and leaves are never collapsed
    assert_eq!(output_nwk, "(A:0.049,B:0.05,C:0.051)root;");
    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_prune_empty_complex_tree() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdgeWithMuts, ()>::new();

    let root = graph.add_node(TestNode(Some("root".to_owned())));
    let internal1 = graph.add_node(TestNode(Some("internal1".to_owned())));
    let internal2 = graph.add_node(TestNode(Some("internal2".to_owned())));
    let internal3 = graph.add_node(TestNode(Some("internal3".to_owned())));
    let a = graph.add_node(TestNode(Some("A".to_owned())));
    let b = graph.add_node(TestNode(Some("B".to_owned())));
    let c = graph.add_node(TestNode(Some("C".to_owned())));
    let d = graph.add_node(TestNode(Some("D".to_owned())));

    // Complex tree structure: root -> internal1 (with muts) -> A (leaf)
    //                              -> internal1 -> internal3 (no muts) -> C,D (leaves)
    //                         root -> internal2 (no muts) -> B (leaf)
    let root_to_internal1 = graph.add_edge(root, internal1, TestEdgeWithMuts::new(Some(0.1), Some(2)))?; // has muts
    let root_to_internal2 = graph.add_edge(root, internal2, TestEdgeWithMuts::new(Some(0.1), Some(0)))?; // no muts, to internal
    let internal1_to_a = graph.add_edge(internal1, a, TestEdgeWithMuts::new(Some(0.1), Some(1)))?; // leaf, has muts
    let internal1_to_internal3 = graph.add_edge(internal1, internal3, TestEdgeWithMuts::new(Some(0.1), Some(0)))?; // no muts, to internal
    let internal2_to_b = graph.add_edge(internal2, b, TestEdgeWithMuts::new(Some(0.1), Some(0)))?; // leaf, no muts (preserved)
    let internal3_to_c = graph.add_edge(internal3, c, TestEdgeWithMuts::new(Some(0.1), Some(1)))?; // leaf, has muts
    let internal3_to_d = graph.add_edge(internal3, d, TestEdgeWithMuts::new(Some(0.1), Some(2)))?; // leaf, has muts

    graph.build()?;

    collapse_short_edges(&mut graph, None, true)?;

    // internal2 and internal3 should be collapsed (empty internal edges), leaves preserved
    // Result: root -> internal1 -> A, C, D and root -> B
    assert_eq!(graph.get_nodes().len(), 6); // root, internal1, A, B, C, D
    assert_eq!(graph.get_edges().len(), 5); // root->internal1, internal1->A, internal1->C, internal1->D, root->B

    Ok(())
  }

  #[test]
  fn test_collapse_short_edges_prune_both_disabled() -> Result<(), Report> {
    let mut graph = Graph::<TestNode, TestEdgeWithMuts, ()>::new();

    let root = graph.add_node(TestNode(Some("root".to_owned())));
    let internal = graph.add_node(TestNode(Some("internal".to_owned())));
    let a = graph.add_node(TestNode(Some("A".to_owned())));

    // Very short edge with no mutations - should be preserved when both pruning options disabled
    let root_to_internal = graph.add_edge(root, internal, TestEdgeWithMuts::new(Some(0.0001), Some(0)))?;
    let internal_to_a = graph.add_edge(internal, a, TestEdgeWithMuts::new(Some(0.1), Some(1)))?;

    graph.build()?;

    collapse_short_edges(&mut graph, None, false)?;

    // Nothing should be collapsed
    assert_eq!(graph.get_nodes().len(), 3); // root, internal, A
    assert_eq!(graph.get_edges().len(), 2); // root->internal, internal->A

    Ok(())
  }
}
