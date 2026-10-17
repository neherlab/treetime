#[cfg(test)]
mod tests {
  use crate::optimize::likelihood::{evaluate_partitions, evaluate_partitions_log_lh_only};
  use crate::pretty_assert_abs_diff_eq;
  use eyre::Report;
  use helpers::{five_point_derivative, setup_skewed_pi};
  use itertools::Itertools;
  use ndarray::Array1;
  use rstest::rstest;

  #[rstest]
  #[case::dense(false)]
  #[case::sparse(true)]
  fn test_branch_objective_matches_tree_likelihood_with_skewed_pi(
    #[case] sparse: bool,
    #[values(1.0, 0.6)] mu: f64,
    #[values(0.02, 0.1, 0.5, 2.0)] t: f64,
  ) -> Result<(), Report> {
    let tree = setup_skewed_pi(sparse, mu)?;
    let (expected, actual): (Vec<f64>, Vec<f64>) = tree
      .graph
      .get_edges()
      .map(|edge| -> Result<(f64, f64), Report> {
        let key = edge.key();
        let contribution = [tree.reconstruction.create_edge_contribution(key)?];
        let objective_0 = evaluate_partitions_log_lh_only(&contribution, tree.branch_lengths[&key])?.value();
        let objective = evaluate_partitions_log_lh_only(&contribution, t)?.value();
        Ok((
          tree.log_lh_with_branch_length(key, t)? - tree.log_lh_0,
          objective - objective_0,
        ))
      })
      .process_results(|changes| changes.unzip())?;
    pretty_assert_abs_diff_eq!(Array1::from_iter(expected), Array1::from_iter(actual), epsilon = 1e-8);
    Ok(())
  }

  #[rstest]
  #[case::dense(false)]
  #[case::sparse(true)]
  fn test_branch_objective_derivative_matches_tree_likelihood_with_skewed_pi(
    #[case] sparse: bool,
    #[values(1.0, 0.6)] mu: f64,
    #[values(0.02, 0.1, 0.5, 2.0)] t: f64,
  ) -> Result<(), Report> {
    let tree = setup_skewed_pi(sparse, mu)?;
    let (expected, actual): (Vec<f64>, Vec<f64>) = tree
      .graph
      .get_edges()
      .map(|edge| -> Result<(f64, f64), Report> {
        let key = edge.key();
        let contribution = [tree.reconstruction.create_edge_contribution(key)?];
        let derivative = evaluate_partitions(&contribution, t)?.derivative;
        let numerical = five_point_derivative(|t| tree.log_lh_with_branch_length(key, t), t)?;
        Ok((numerical, derivative))
      })
      .process_results(|derivatives| derivatives.unzip())?;
    pretty_assert_abs_diff_eq!(Array1::from_iter(expected), Array1::from_iter(actual), epsilon = 1e-9);
    Ok(())
  }

  #[rstest]
  #[case::dense(false)]
  #[case::sparse(true)]
  fn test_branch_objective_second_derivative_matches_derivative_change_with_skewed_pi(
    #[case] sparse: bool,
    #[values(1.0, 0.6)] mu: f64,
    #[values(0.02, 0.1, 0.5, 2.0)] t: f64,
  ) -> Result<(), Report> {
    let tree = setup_skewed_pi(sparse, mu)?;
    let (expected, actual): (Vec<f64>, Vec<f64>) = tree
      .graph
      .get_edges()
      .map(|edge| -> Result<(f64, f64), Report> {
        let contribution = [tree.reconstruction.create_edge_contribution(edge.key())?];
        let second_derivative = evaluate_partitions(&contribution, t)?.second_derivative;
        let numerical = five_point_derivative(|t| Ok(evaluate_partitions(&contribution, t)?.derivative), t)?;
        Ok((numerical, second_derivative))
      })
      .process_results(|derivatives| derivatives.unzip())?;
    pretty_assert_abs_diff_eq!(Array1::from_iter(expected), Array1::from_iter(actual), epsilon = 1e-8);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::gtr::gtr::{GTR, Mu};
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::dense::partition::PartitionMarginalDense;
    use crate::partition::marginal::reconstruction::{
      DenseReconstruction, MarginalReconstruction, SparseReconstruction,
    };
    use crate::test_utils::leaf_seq_inputs;
    use eyre::Report;
    use indoc::indoc;
    use ndarray::array;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_io::fasta::fasta_read;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;

    pub(super) struct SkewedPiTree {
      pub(super) graph: Graph,
      pub(super) branch_lengths: BTreeMap<GraphEdgeKey, f64>,
      pub(super) reconstruction: MarginalReconstruction,
      pub(super) log_lh_0: f64,
    }

    impl SkewedPiTree {
      pub(super) fn log_lh_with_branch_length(&self, key: GraphEdgeKey, t: f64) -> Result<f64, Report> {
        let mut changed = self.branch_lengths.clone();
        changed.insert(key, t);
        let (_, log_lh) = self.reconstruction.clone().marginal_update(&self.graph, &changed)?;
        Ok(log_lh.value())
      }
    }

    pub(super) fn five_point_derivative(f: impl Fn(f64) -> Result<f64, Report>, t: f64) -> Result<f64, Report> {
      let h = 1e-3 * t;
      Ok((f(t - 2.0 * h)? - 8.0 * f(t - h)? + 8.0 * f(t + h)? - f(t + 2.0 * h)?) / (12.0 * h))
    }

    pub(super) fn setup_skewed_pi(sparse: bool, mu: f64) -> Result<SkewedPiTree, Report> {
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
      let gtr = GTR::builder()
        .n_states(alphabet.n_canonical())
        .mu(Mu::Normalized(mu))
        .pi(array![0.7, 0.2, 0.07, 0.03])
        .build()?;

      let inputs = leaf_seq_inputs(&graph, &names, aln);
      let reconstruction = if sparse {
        let (partition, node_states) =
          create_fitch_partition(&graph, alphabet, inputs)?.into_marginal_sparse(&graph)?;
        MarginalReconstruction::Sparse(SparseReconstruction::seeded(partition, gtr, node_states))
      } else {
        let partition = PartitionMarginalDense::new(alphabet, &graph, &inputs)?;
        MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, gtr))
      };
      let (reconstruction, log_lh_0) = reconstruction.marginal_update(&graph, &branch_lengths)?;

      Ok(SkewedPiTree {
        graph,
        branch_lengths,
        reconstruction,
        log_lh_0: log_lh_0.value(),
      })
    }
  }
}
