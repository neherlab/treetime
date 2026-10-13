#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::gtr::GTR;
  use crate::optimize::likelihood::evaluate_mixed_log_lh_only;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
  use crate::pretty_assert_abs_diff_eq;
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  #[rstest]
  #[case::dense(false)]
  #[case::sparse(true)]
  fn test_branch_objective_matches_tree_likelihood_with_skewed_pi(
    #[case] sparse: bool,
    #[values(0.02, 0.1, 0.5, 2.0)] t: f64,
  ) -> Result<(), Report> {
    let parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let branch_lengths: BTreeMap<GraphEdgeKey, f64> = parsed
      .branch_lengths
      .iter()
      .map(|(&k, v)| (k, v.unwrap_or(0.0)))
      .collect();
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = fasta_read(
      indoc! {r#"
      >A
      AAAAAAAACATAAAGA
      >B
      AAAAAACAAATAAAGA
      >C
      AAAATAAAACAAAGAA
      >D
      AAAATAAGACAAATAA
      "#}
      .as_bytes(),
      &alphabet,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
    let mut gtr = GTR::builder()
      .n_states(alphabet.n_canonical())
      .mu(1.0)
      .pi(array![0.7, 0.2, 0.07, 0.03])
      .build()?;
    // TODO: the branch objective omits mu (kb/issues/M-branch-objective-omits-substitution-rate.md), so the test
    // needs mu = 1
    gtr.mu = 1.0;

    let inputs = leaf_seq_inputs(&graph, &names, aln);
    let reconstruction = if sparse {
      let (partition, node_states) = create_fitch_partition(&graph, alphabet, inputs)?.into_marginal_sparse(&graph)?;
      MarginalReconstruction::Sparse(SparseReconstruction::seeded(partition, gtr, node_states))
    } else {
      let partition = PartitionMarginalDense::new(alphabet, &graph, &inputs)?;
      MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, gtr))
    };
    let (reconstruction, log_lh_0) = reconstruction.marginal_update(&graph, &branch_lengths)?;

    let (expected, actual): (Vec<f64>, Vec<f64>) = graph
      .get_edges()
      .map(|edge| -> Result<(f64, f64), Report> {
        let key = edge.key();
        let contribution = [reconstruction.create_edge_contribution(key)?];
        let objective_0 = evaluate_mixed_log_lh_only(&contribution, branch_lengths[&key])?.value();
        let objective = evaluate_mixed_log_lh_only(&contribution, t)?.value();
        let mut changed = branch_lengths.clone();
        changed.insert(key, t);
        let (_, log_lh) = reconstruction.clone().marginal_update(&graph, &changed)?;
        Ok((log_lh.value() - log_lh_0.value(), objective - objective_0))
      })
      .process_results(|changes| changes.unzip())?;
    pretty_assert_abs_diff_eq!(
      Array1::from_iter(expected),
      Array1::from_iter(actual),
      epsilon = 1e-8
    );
    Ok(())
  }
}
