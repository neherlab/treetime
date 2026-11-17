#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::gamma::GammaRates;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::partition_model::PartitionModel;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, SparseReconstruction};
  use crate::pretty_assert_abs_diff_eq;
  use crate::test_utils::{find_node_key_by_name, leaf_seq_inputs};
  use eyre::Report;
  use helpers::expected_leaf_posterior;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, AsciiChar, Seq};

  const TREE: &str = "(A:0.1,B:0.2)root;";

  #[rustfmt::skip]
  #[rstest]
  #[case::missing(  0, b'N')]
  #[case::ambiguous(1, b'R')]
  #[trace]
  fn test_mixture_leaf_posteriors_weight_the_category_posteriors(
    #[case] pos: usize,
    #[case] observed: u8,
  ) -> Result<(), Report> {
    let parsed = nwk_read(TREE.as_bytes())?;
    let graph = parsed.graph;
    let branch_lengths = branch_lengths_or_zero(&parsed.branch_lengths);
    let records = [("A", "AAC"), ("B", "NRC")]
      .into_iter()
      .map(|(name, seq)| {
        Ok(AlignmentRecord {
          name: name.to_owned(),
          seq: Seq::try_from_str(seq)?,
        })
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let inputs = leaf_seq_inputs(&graph, &parsed.names, records);
    let gamma = GammaRates::new(4, 0.5)?;
    let model = PartitionModel::new(jc69(JC69Params::default())?, gamma.clone());
    let leaf = find_node_key_by_name(&graph, &parsed.names, "B").expect("leaf B exists");
    let (root, edge_key) = graph.node_parent(leaf)?.expect("leaf B has a parent");
    let expected = expected_leaf_posterior(&gamma, observed, 0.3);

    let dense = DenseReconstruction::seeded(PartitionMarginalDense::new(Alphabet::default(), &graph, &inputs)?, model.clone())
      .marginal_update(&graph, &branch_lengths)?
      .0;
    let dense_actual = dense.node_states[&leaf].profile.dis.row(pos).to_owned();

    let (partition, nodes) = create_fitch_partition(&graph, Alphabet::default(), inputs)?.into_marginal_sparse(&graph)?;
    let sparse = SparseReconstruction::seeded(partition, model, nodes)
      .marginal_update(&graph, &branch_lengths)?
      .0;
    let sparse_actual = sparse
      .leaf_downs()
      .down(edge_key, root)
      .posterior(AsciiChar::from_byte_unchecked(observed), pos, &Alphabet::default())
      .expect("an unresolved leaf state has a posterior");

    pretty_assert_abs_diff_eq!(expected, dense_actual, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(expected, sparse_actual, epsilon = 1e-12);
    Ok(())
  }

  mod helpers {
    use crate::gtr::gamma::GammaRates;
    use ndarray::{Array1, array};

    pub(super) fn expected_leaf_posterior(gamma: &GammaRates, observed: u8, path_length: f64) -> Array1<f64> {
      let mask = match observed {
        b'R' => array![1.0, 0.0, 1.0, 0.0],
        _ => array![1.0, 1.0, 1.0, 1.0],
      };
      let category = |rate: f64| -> Array1<f64> {
        let decay = (-4.0 / 3.0 * rate * path_length).exp();
        let same = 0.25 + 0.75 * decay;
        let other = 0.25 - 0.25 * decay;
        array![same, other, other, other] * &mask
      };
      let total: Array1<f64> = gamma
        .rates()
        .iter()
        .map(|&rate| category(rate))
        .fold(Array1::zeros(4), |acc, p| acc + p);
      let norm = total.sum();
      total / norm
    }
  }
}
