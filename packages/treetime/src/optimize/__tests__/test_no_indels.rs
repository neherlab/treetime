#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::args::{BranchOptMethod, InitialGuessMode};
  use crate::optimize::iteration::apply_initial_guess_mode;
  use crate::optimize::optimize_unified::run_optimize_mixed_inner;
  use crate::optimize::run_loop::run_optimize_loop;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  const TREE_ZERO_BL: &str = "((A:0.0,B:0.0)AB:0.0,C:0.0)root:0.0;";

  fn setup_dense_with_marginal(
    newick: &str,
  ) -> Result<(GraphAncestral, Vec<Arc<RwLock<PartitionMarginalDense>>>), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let aln = read_many_fasta_str(">A\nACGT\n>B\nACGT\n>C\nACGT\n", &alphabet)?;
    let gtr = jc69(JC69Params::default())?;
    let partition = PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    let partitions = vec![Arc::new(RwLock::new(partition))];
    initialize_marginal(&graph, &partitions, &aln)?;
    update_marginal(&graph, &partitions)?;
    Ok((graph, partitions))
  }

  fn inject_indel_on_first_edge(
    graph: &GraphAncestral,
    partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  ) -> Result<(), Report> {
    let first_edge_key = graph.get_edges()[0].read_arc().key();
    partitions[0].write_arc().edges.get_mut(&first_edge_key).unwrap().indels =
      vec![InDel::del((0, 2), Seq::try_from_str("AC")?)];
    Ok(())
  }

  #[test]
  fn test_no_indels_drops_indel_contribution_from_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph_with: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_with, sparse_with, mixed_with) = setup_partitions(&graph_with, &aln)?;

    let mut graph_without: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_without, sparse_without, mixed_without) = setup_partitions(&graph_without, &aln)?;

    let first_edge_key = graph_with.get_edges()[0].read_arc().key();
    graph_with.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.1));
    sparse_with[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];

    let first_edge_key_without = graph_without.get_edges()[0].read_arc().key();
    graph_without.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.1));
    sparse_without[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key_without)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];

    let result_with = run_optimize_loop(
      &mut graph_with,
      &sparse_with,
      &dense_with,
      &mixed_with,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
    )?;

    let result_without = run_optimize_loop(
      &mut graph_without,
      &sparse_without,
      &dense_without,
      &mixed_without,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      true,
    )?;

    assert!(
      result_without.lh_history[0] > result_with.lh_history[0],
      "no_indels=true should produce higher likelihood than with indels: {} vs {}",
      result_without.lh_history[0],
      result_with.lh_history[0]
    );
    Ok(())
  }

  #[test]
  fn test_no_indels_optimizer_ignores_indel_counts() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (_, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let first_edge_key = graph.get_edges()[0].read_arc().key();
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.05));
    sparse_partitions[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];

    update_marginal(&graph, &sparse_partitions)?;

    let bl_before: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect();

    run_optimize_mixed_inner(&graph, &mixed_partitions, BranchOptMethod::BrentSqrt, 0.0, true)?;

    let bl_after: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect();

    assert_ne!(bl_before, bl_after, "Optimizer should modify branch lengths");
    Ok(())
  }

  #[test]
  fn test_no_indels_matches_no_indel_data() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph_no_flag: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_nf, sparse_nf, mixed_nf) = setup_partitions(&graph_no_flag, &aln)?;

    let mut graph_flag: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_f, sparse_f, mixed_f) = setup_partitions(&graph_flag, &aln)?;

    let result_no_flag = run_optimize_loop(
      &mut graph_no_flag,
      &sparse_nf,
      &dense_nf,
      &mixed_nf,
      3,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
    )?;

    let result_flag = run_optimize_loop(
      &mut graph_flag,
      &sparse_f,
      &dense_f,
      &mixed_f,
      3,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      true,
    )?;

    assert_eq!(
      result_no_flag.lh_history.len(),
      result_flag.lh_history.len(),
      "Same number of iterations expected"
    );
    for (i, (lh_nf, lh_f)) in result_no_flag
      .lh_history
      .iter()
      .zip(&result_flag.lh_history)
      .enumerate()
    {
      assert_abs_diff_eq!(lh_nf, lh_f, epsilon = 1e-10);
    }
    Ok(())
  }

  #[test]
  fn test_no_indels_initial_guess_never_accepts_zero_bl_with_indels() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    let result = apply_initial_guess_mode(&graph, &partitions, InitialGuessMode::Never, true);
    assert!(
      result.is_ok(),
      "no_indels=true should accept zero-BL indel edges in Never mode, got: {result:?}"
    );
    Ok(())
  }
}
