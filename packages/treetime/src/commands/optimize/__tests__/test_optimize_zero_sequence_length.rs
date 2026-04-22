#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::optimize::args::BranchOptMethod;
  use crate::commands::optimize::optimize_unified::{
    initial_guess_mixed, run_optimize_mixed, run_optimize_mixed_with_indel_rate,
  };
  use crate::commands::optimize::run::collect_optimize_partitions;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

  fn zero_length_partitions(graph: &GraphAncestral) -> crate::commands::optimize::partition_ops::PartitionOptimizeVec {
    let dense = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default()).unwrap(),
      alphabet: Alphabet::new(AlphabetName::Nuc).unwrap(),
      length: 0,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];
    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    collect_optimize_partitions(&dense, &sparse)
  }

  #[test]
  fn test_optimize_zero_sequence_length_run_optimize_error() {
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;").unwrap();
    let partitions = zero_length_partitions(&graph);
    let result = run_optimize_mixed(&graph, &partitions, BranchOptMethod::Newton);
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot optimize branch lengths"
    );
  }

  #[test]
  fn test_optimize_zero_sequence_length_initial_guess_error() {
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;").unwrap();
    let partitions = zero_length_partitions(&graph);
    let result = initial_guess_mixed(&graph, &partitions, true);
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot compute initial guess"
    );
  }

  #[test]
  fn test_optimize_zero_sequence_length_run_optimize_with_fixed_rate_error() {
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;").unwrap();
    let partitions = zero_length_partitions(&graph);
    let result = run_optimize_mixed_with_indel_rate(&graph, &partitions, BranchOptMethod::Newton, 1.0);
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot optimize branch lengths"
    );
  }
}
