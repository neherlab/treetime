use crate::alphabet::alphabet::Alphabet;
use crate::gtr::get_gtr::{GtrModelName, get_gtr_by_name, log_gtr};
use crate::gtr::gtr::GTR;
use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::marginal_sparse::PartitionMarginalSparse;
use crate::payload::ancestral::GraphAncestral;
use crate::prune::prune::prune_nodes;
use eyre::Report;
use parking_lot::RwLock;
use serde::Serialize;
use std::collections::BTreeSet;
use std::sync::Arc;
use treetime_graph::assign_node_names::assign_node_names;
use treetime_io::fasta::FastaRecord;

pub struct PruneParams {
  pub prune_short: Option<f64>,
  pub prune_empty: bool,
  pub merge_shared_mutations: bool,
  pub node_names: BTreeSet<String>,
}

pub struct PruneInput {
  pub graph: GraphAncestral,
  pub alphabet: Alphabet,
  pub sequences: Option<Vec<FastaRecord>>,
}

#[derive(Debug, Serialize)]
pub struct PruneOutput {
  #[serde(skip)]
  pub graph: GraphAncestral,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>>,
}

pub fn run(params: &PruneParams, mut input: PruneInput) -> Result<PruneOutput, Report> {
  let needs_sequences = params.prune_empty || params.merge_shared_mutations;
  let partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = if needs_sequences {
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
    )?;
    match created.partition {
      MarginalPartition::Sparse(p) => vec![Arc::new(RwLock::new(p))],
      MarginalPartition::Dense(_) => {
        let gtr = get_gtr_by_name(GtrModelName::JC69)?;
        log_gtr(&gtr, GtrModelName::JC69);
        let fitch =
          crate::ancestral::fitch::create_fitch_partition(&input.graph, 0, input.alphabet.clone(), sequences)?;
        let partition = fitch.into_marginal_sparse(gtr, &input.graph)?;
        vec![Arc::new(RwLock::new(partition))]
      },
    }
  } else {
    vec![]
  };

  prune_nodes(
    &mut input.graph,
    &partitions,
    params.prune_short,
    params.prune_empty,
    &params.node_names,
  )?;

  if params.merge_shared_mutations {
    merge_shared_mutation_branches(&mut input.graph, &partitions)?;
    input.graph.build()?;
    assign_node_names(&input.graph)?;
  }

  let gtr = (!partitions.is_empty()).then(|| partitions[0].read_arc().gtr.clone());

  Ok(PruneOutput {
    graph: input.graph,
    gtr,
    partitions,
  })
}
