use crate::alphabet::alphabet::Alphabet;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr};
use crate::gtr::gtr::GTR;
use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::prune::prune::prune_nodes;
use eyre::Report;
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::assign_node_names::assign_node_names;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaRecord;

pub struct PruneParams {
  pub prune_short: Option<f64>,
  pub prune_empty: bool,
  pub merge_shared_mutations: bool,
  pub node_names: BTreeSet<String>,
}

pub struct PruneInput {
  pub graph: Graph,
  pub alphabet: Alphabet,
  pub sequences: Option<Vec<FastaRecord>>,
  /// Raw per-edge branch lengths captured from the Newick parse, keyed by edge id. The collapse and
  /// merge producers update it in place across the topology edits; it exits as
  /// `PruneOutput.branch_lengths`.
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

#[derive(Debug, Serialize)]
pub struct PruneOutput {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub partitions: Vec<PartitionMarginalSparse>,
  #[serde(skip)]
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  #[serde(skip)]
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
}

pub fn run(
  params: &PruneParams,
  mut input: PruneInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<PruneOutput, Report> {
  // Entry maps propagated through the whole prune pipeline: every downstream name and branch length
  // read comes from these maps, not the payload. `names` comes from the parse; `assign_node_names`
  // after a merge refreshes it to the post-topology labels; the collapse and merge producers maintain
  // `branch_lengths` in place across the topology edits, so both maps exit reflecting the final pruned
  // tree.
  let mut names = names.clone();
  let mut branch_lengths = std::mem::take(&mut input.branch_lengths);

  let needs_sequences = params.prune_empty || params.merge_shared_mutations;
  let mut partitions: Vec<PartitionMarginalSparse> = if needs_sequences {
    let sequences = input
      .sequences
      .as_ref()
      .ok_or_else(|| eyre::eyre!("Sequences required for --prune-empty or --merge-shared-mutations"))?;
    let created = create_marginal_partition(
      &input.graph,
      0,
      input.alphabet.clone(),
      sequences,
      GtrModelName::JC69,
      None,
      &branch_lengths,
      &names,
    )?;
    match created.partition {
      MarginalPartition::Sparse(p) => vec![p],
      MarginalPartition::Dense(_) => {
        let gtr = get_gtr_by_name(GtrModelName::JC69)?;
        log_gtr(&gtr, GtrModelName::JC69);
        let fitch =
          crate::ancestral::fitch::create_fitch_partition(&input.graph, 0, input.alphabet.clone(), sequences, &names)?;
        let partition = fitch.into_marginal_sparse(gtr, &input.graph)?;
        vec![partition]
      },
    }
  } else {
    vec![]
  };

  prune_nodes(
    &mut input.graph,
    &mut partitions,
    params.prune_short,
    params.prune_empty,
    &params.node_names,
    &names,
    &mut branch_lengths,
  )?;

  if params.merge_shared_mutations {
    // The merge producer updates the branch-length map in place; refresh `names` from the post-merge
    // topology so downstream name readers get the reassigned labels.
    merge_shared_mutation_branches(&mut input.graph, &mut partitions, &mut branch_lengths)?;
    input.graph.build()?;
    names = assign_node_names(names, &input.graph)?;
  }

  let gtr = (!partitions.is_empty()).then(|| partitions[0].gtr.clone());

  Ok(PruneOutput {
    graph: input.graph,
    gtr,
    partitions,
    names,
    branch_lengths,
  })
}
