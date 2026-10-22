use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeAminoAcids, TreeSequences};
use app_output::augur_node_data_ancestral::{AncestralNodeSequences, write_augur_node_data_ancestral};
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::tree_output::{write_partition_tree_outputs, write_run_tree_outputs};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime::alphabet::alphabet::Alphabet;
use treetime::ancestral::aa::AaNodeData;
use treetime::progress::LogSink;
use treetime::progress_info;
use treetime::seq::mutation::{Mutation, PositionBase};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node_names::NodeNames;
use treetime_primitives::Seq;
use treetime_utils::make_internal_error;
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub(crate) fn write_ancestral_trees(
  trees: &AncestralTrees<'_>,
  partitions: Vec<AncestralPartitionTrees<'_>>,
  resolved: &ResolvedOutputs,
  command: CommandKind,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let numbers = partitions.iter().map(|partition| partition.number).collect_vec();
  if !numbers.iter().copied().eq(resolved.partitions()) {
    return make_internal_error!(
      "Ancestral outputs are planned for {} partition(s), but the reconstruction returned partitions {numbers:?}",
      resolved.partition_count
    );
  }
  let run_graph = AnnotatedGraph {
    graph: trees.graph,
    names: trees.names,
    divergence_branch_lengths: trees.branch_lengths,
    time_branch_lengths: None,
    divergence: Divergence::CumulativeBranchLength,
    sequences: None,
    dates: None,
    traits: None,
  };
  write_run_tree_outputs(&run_graph, resolved, command, log)?;
  for partition in partitions {
    let annotated = AnnotatedGraph {
      sequences: Some(TreeSequences {
        alphabet: partition.alphabet,
        root_sequence: &partition.maps.root_sequence,
        edge_mutations: &partition.maps.edge_mutations,
        amino_acids: trees
          .amino_acids
          .filter(|_| partition.number == 1)
          .map(|(node_data, cdses)| TreeAminoAcids { node_data, cdses }),
        comment_positions: trees.comment_positions,
      }),
      ..run_graph
    };
    let tree = write_partition_tree_outputs(&annotated, partition.number, resolved, command, log)?;
    if let (Some(tree), Some(path), Some(node_sequences)) = (
      &tree,
      resolved.partition_path(OutputSelection::NodeData, partition.number),
      partition.node_sequences,
    ) {
      write_augur_node_data_ancestral(tree, node_sequences, path)?;
      progress_info!(log, "Wrote augur node data JSON to {}", path.display());
    }
  }
  Ok(())
}

pub(crate) struct AncestralTrees<'a> {
  pub(crate) graph: &'a Graph,
  pub(crate) names: &'a NodeNames,
  pub(crate) branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  pub(crate) amino_acids: Option<&'a (AaNodeData, BTreeMap<String, AugurNodeDataJsonAnnotationEntry>)>,
  pub(crate) comment_positions: PositionBase,
}

pub(crate) struct AncestralPartitionTrees<'a> {
  pub(crate) number: usize,
  pub(crate) alphabet: &'a Alphabet,
  pub(crate) maps: AncestralOutputMaps,
  pub(crate) node_sequences: Option<AncestralNodeSequences<'a>>,
}

pub(crate) struct AncestralOutputMaps {
  pub(crate) root_sequence: Seq,
  pub(crate) edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
}
