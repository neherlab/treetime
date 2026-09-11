#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::{initial_guess_mixed, run_optimize_mixed, run_optimize_mixed_with_indel_rate};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::optimize_partition_view;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::payload::ancestral::GraphAncestral;

  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_utils::assert_error;

  fn zero_length_partitions(_graph: &GraphAncestral) -> (Vec<PartitionMarginalDense>, Vec<PartitionMarginalSparse>) {
    let dense = vec![PartitionMarginalDense::new(
      0,
      jc69(JC69Params::default()).unwrap(),
      Alphabet::new(AlphabetName::Nuc).unwrap(),
      0,
    )];
    let sparse: Vec<PartitionMarginalSparse> = vec![];
    (dense, sparse)
  }

  #[test]
  fn test_optimize_zero_sequence_length_run_optimize_error() {
    let NwkParse {
      graph,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;").unwrap();
    let graph: GraphAncestral = graph;
    let (dense, sparse) = zero_length_partitions(&graph);
    let partitions = optimize_partition_view(&dense, &sparse);
    let result = run_optimize_mixed(&graph, &partitions, BranchOptMethod::Newton, &mut branch_lengths);
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot optimize branch lengths"
    );
  }

  #[test]
  fn test_optimize_zero_sequence_length_initial_guess_error() {
    let NwkParse {
      graph,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;").unwrap();
    let graph: GraphAncestral = graph;
    let (dense, sparse) = zero_length_partitions(&graph);
    let partitions = optimize_partition_view(&dense, &sparse);
    let result = initial_guess_mixed(&graph, &partitions, true, false, &mut branch_lengths);
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot compute initial guess"
    );
  }

  #[test]
  fn test_optimize_zero_sequence_length_run_optimize_with_fixed_rate_error() {
    let NwkParse {
      graph,
      mut branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;").unwrap();
    let graph: GraphAncestral = graph;
    let (dense, sparse) = zero_length_partitions(&graph);
    let partitions = optimize_partition_view(&dense, &sparse);
    let result =
      run_optimize_mixed_with_indel_rate(&graph, &partitions, BranchOptMethod::Newton, 1.0, &mut branch_lengths);
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot optimize branch lengths"
    );
  }
}
