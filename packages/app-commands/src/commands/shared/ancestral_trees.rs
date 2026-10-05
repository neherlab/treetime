use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeAminoAcids, TreeSequences};
use app_output::augur_node_data_ancestral::{AncestralNodeSequences, write_augur_node_data_ancestral};
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
use std::collections::BTreeMap;
use treetime::ancestral::aa::AaNodeData;
use treetime::progress::LogSink;
use treetime::progress_info;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub(crate) struct AncestralTrees<'a> {
  pub(crate) graph: &'a Graph,
  pub(crate) names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  pub(crate) branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  pub(crate) maps: &'a AncestralOutputMaps,
  pub(crate) amino_acids: Option<&'a (AaNodeData, BTreeMap<String, AugurNodeDataJsonAnnotationEntry>)>,
}

pub(crate) struct AncestralOutputMaps {
  pub(crate) root_sequence: Seq,
  pub(crate) edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}

pub(crate) fn write_ancestral_trees(
  trees: &AncestralTrees<'_>,
  node_sequences: Option<AncestralNodeSequences<'_>>,
  resolved: &ResolvedOutputs,
  command: CommandKind,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let annotated = AnnotatedGraph {
    graph: trees.graph,
    names: trees.names,
    divergence_branch_lengths: trees.branch_lengths,
    time_branch_lengths: None,
    divergence: Divergence::CumulativeBranchLength,
    sequences: Some(TreeSequences {
      root_sequence: &trees.maps.root_sequence,
      edge_mutations: &trees.maps.edge_mutations,
      mutation_counts: None,
      amino_acids: trees
        .amino_acids
        .map(|(node_data, cdses)| TreeAminoAcids { node_data, cdses }),
    }),
    dates: None,
    traits: None,
  };
  let tree = tree_view_for_outputs(&annotated, resolved);
  if let (Ok(Some(tree)), Some(path), Some(node_sequences)) =
    (&tree, resolved.path(OutputSelection::AugurNodeData), node_sequences)
  {
    write_augur_node_data_ancestral(tree, node_sequences, path)?;
    progress_info!(log, "Wrote augur node data JSON to {}", path.display());
  }
  write_graph_outputs(&annotated, &resolved.tree_outputs)?;
  let Some(tree) = tree? else {
    return Ok(());
  };
  write_tree_outputs(&tree, &resolved.tree_outputs, command, log)
}
