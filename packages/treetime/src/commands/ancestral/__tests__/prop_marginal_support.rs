#![cfg(test)]

use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::commands::ancestral::__tests__::prop_generators::input::MarginalTestInput;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::payload::ancestral::GraphAncestral;
use eyre::Report;
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_io::nwk::nwk_read_str;

pub fn run_dense_marginal(
  input: &MarginalTestInput,
) -> Result<(f64, [Arc<RwLock<PartitionMarginalDense>>; 1]), Report> {
  let graph: GraphAncestral = nwk_read_str(&input.newick)?;
  let alphabet = Alphabet::new(AlphabetName::Nuc, true)?;
  let length = get_common_length(&input.alignment)?;

  let partitions = [Arc::new(RwLock::new(PartitionMarginalDense {
    index: 0,
    gtr: input.gtr.clone(),
    alphabet,
    length,
    nodes: btreemap! {},
    edges: btreemap! {},
  }))];

  let log_lh = initialize_marginal(&graph, &partitions, &input.alignment)?;
  Ok((log_lh, partitions))
}

pub fn run_sparse_marginal(
  input: &MarginalTestInput,
) -> Result<(f64, [Arc<RwLock<PartitionMarginalSparse>>; 1]), Report> {
  let graph: GraphAncestral = nwk_read_str(&input.newick)?;
  let alphabet = Alphabet::default();
  let length = get_common_length(&input.alignment)?;

  let partitions = [Arc::new(RwLock::new(PartitionMarginalSparse {
    index: 0,
    gtr: input.gtr.clone(),
    alphabet,
    length,
    nodes: btreemap! {},
    edges: btreemap! {},
  }))];

  compress_sequences(&graph, &partitions, &input.alignment)?;
  let log_lh = update_marginal(&graph, &partitions)?;
  Ok((log_lh, partitions))
}
