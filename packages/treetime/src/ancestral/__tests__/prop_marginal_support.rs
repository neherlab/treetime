#[cfg(test)]
pub(super) mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::__tests__::prop_generators::input::MarginalTestInput;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;

  use treetime_graph::graph::Graph;

  use treetime_io::nwk::nwk_read_str;

  pub(crate) fn run_dense_marginal(input: &MarginalTestInput) -> Result<(f64, DenseReconstruction), Report> {
    let nwk_parsed = nwk_read_str(&input.newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let length = get_common_length(&input.alignment)?;

    let partition = PartitionMarginalDense::new(0, alphabet, length);
    let node_states = partition.attach_sequences(&graph, &node_seq_inputs(&graph, &names, input.alignment.clone()))?;
    let recon = DenseReconstruction::seeded(partition, input.gtr.clone(), node_states);
    let (recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }

  pub(crate) fn run_sparse_marginal(input: &MarginalTestInput) -> Result<(f64, SparseReconstruction), Report> {
    let nwk_parsed = nwk_read_str(&input.newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let alphabet = Alphabet::default();
    let _ = get_common_length(&input.alignment)?;

    let fitch = create_fitch_partition(
      &graph,
      0,
      alphabet,
      &node_seq_inputs(&graph, &names, input.alignment.clone()),
    )?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, input.gtr.clone(), node_states);
    let (recon, log_lh) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let log_lh = log_lh.value();
    Ok((log_lh, recon))
  }
}
