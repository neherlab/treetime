#[cfg(test)]
mod tests {
  use crate::gtr::gamma::GammaRates;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::partition_model::PartitionModel;
  use crate::partition::marginal::mixture::__tests__::fixtures::tests::helpers::{
    GAP_FREE, LONG_BRANCHES, MISSING_DATA, fixture, scaled_fixture,
  };
  use crate::partition::marginal::mixture::contribution::{dense_mixture_contribution, sparse_mixture_contribution};
  use crate::partition::marginal::mixture::dense::dense_mixture_update;
  use crate::partition::marginal::mixture::sparse::sparse_mixture_update;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::pretty_assert_abs_diff_eq;
  use crate::test_utils::make_nonuniform_gtr;
  use eyre::Report;
  use helpers::{
    dense_profiles, dense_sequences, edge_indels, sparse_estimates, sparse_posteriors, sparse_sites, stacked,
  };
  use itertools::Itertools;
  use ndarray::array;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::gap_free(    &GAP_FREE)]
  #[case::missing_data(&MISSING_DATA)]
  #[trace]
  fn test_dense_single_category_mixture_matches_the_plain_pass(#[case] alignment: &[(&str, &str)]) -> Result<(), Report> {
    let fixture = fixture(alignment)?;
    let gtr = make_nonuniform_gtr()?;
    let expected = fixture.dense.marginal_update(&gtr, &fixture.graph, &fixture.branch_lengths, &())?;

    let model = PartitionModel::from(gtr.clone());
    let (actual, state) = dense_mixture_update(&fixture.dense, &model, &fixture.graph, &fixture.branch_lengths)?;

    let edges = fixture.graph.get_edges().map(|edge| edge.key()).collect_vec();
    let plain_contributions = edges
      .iter()
      .map(|&edge_key| {
        fixture
          .dense
          .create_edge_contribution(&gtr, &expected.edges.backward, &expected.edges.forward, edge_key)
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let mixture_contributions = edges
      .iter()
      .map(|&edge_key| dense_mixture_contribution(&fixture.dense, &state, &array![1.0], edge_key))
      .collect::<Result<Vec<_>, Report>>()?;

    pretty_assert_abs_diff_eq!(expected.log_lh.value(), actual.log_lh.value(), epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(dense_profiles(&expected.node_states)?, dense_profiles(&actual.node_states)?, epsilon = 1e-12);
    assert_eq!(dense_sequences(&expected.node_states), dense_sequences(&actual.node_states));
    pretty_assert_abs_diff_eq!(stacked(&plain_contributions)?, stacked(&mixture_contributions)?, epsilon = 1e-12);
    assert_eq!(edge_indels(&expected.edges.estimates, &edges), edge_indels(&actual.edges.estimates, &edges));
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::gap_free(    &GAP_FREE)]
  #[case::missing_data(&MISSING_DATA)]
  #[trace]
  fn test_sparse_single_category_mixture_matches_the_plain_pass(#[case] alignment: &[(&str, &str)]) -> Result<(), Report> {
    let fixture = fixture(alignment)?;
    let gtr = make_nonuniform_gtr()?;
    let expected =
      fixture
        .sparse
        .marginal_update(&gtr, &fixture.graph, &fixture.branch_lengths, &fixture.sparse_nodes)?;

    let model = PartitionModel::from(gtr.clone());
    let (actual, state) = sparse_mixture_update(
      &fixture.sparse,
      &model,
      &fixture.graph,
      &fixture.branch_lengths,
      &fixture.sparse_nodes,
    )?;

    let edges = fixture.graph.get_edges().map(|edge| edge.key()).collect_vec();
    let plain_sites = edges
      .iter()
      .map(|&edge_key| {
        let contribution =
          fixture
            .sparse
            .create_edge_contribution(&gtr, &expected.edges.backward, &expected.edges.forward, edge_key)?;
        Ok(sparse_sites(&contribution))
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let mixture_sites = edges
      .iter()
      .map(|&edge_key| {
        let contribution =
          sparse_mixture_contribution(&fixture.sparse, &actual.node_states, &state, &array![1.0], edge_key)?;
        Ok(sparse_sites(&contribution))
      })
      .collect::<Result<Vec<_>, Report>>()?;
    let length = fixture.sparse.length;

    pretty_assert_abs_diff_eq!(expected.log_lh.value(), actual.log_lh.value(), epsilon = 1e-10);
    pretty_assert_abs_diff_eq!(sparse_posteriors(&expected.node_states, length)?, sparse_posteriors(&actual.node_states, length)?, epsilon = 1e-12);
    assert_eq!(sparse_estimates(&expected.edges.estimates, &edges), sparse_estimates(&actual.edges.estimates, &edges));
    assert_eq!(plain_sites, mixture_sites);
    Ok(())
  }

  #[test]
  fn test_dense_and_sparse_gamma_mixtures_agree() -> Result<(), Report> {
    let fixture = scaled_fixture(&GAP_FREE, LONG_BRANCHES)?;
    let model = PartitionModel::new(make_nonuniform_gtr()?, GammaRates::new(4, 5.0)?);
    let (dense, _) = dense_mixture_update(&fixture.dense, &model, &fixture.graph, &fixture.branch_lengths)?;
    let (sparse, _) = sparse_mixture_update(
      &fixture.sparse,
      &model,
      &fixture.graph,
      &fixture.branch_lengths,
      &fixture.sparse_nodes,
    )?;

    pretty_assert_abs_diff_eq!(dense.log_lh.value(), sparse.log_lh.value(), epsilon = 1e-8);
    pretty_assert_abs_diff_eq!(
      dense_profiles(&dense.node_states)?,
      sparse_posteriors(&sparse.node_states, fixture.sparse.length)?,
      epsilon = 1e-8
    );
    Ok(())
  }

  #[test]
  fn test_gamma_mixture_has_no_zero_branch_certificate() -> Result<(), Report> {
    let fixture = fixture(&GAP_FREE)?;
    let jc69 = jc69(JC69Params::default())?;
    let model = PartitionModel::new(jc69.clone(), GammaRates::new(4, 0.5)?);
    let (_, state) = dense_mixture_update(&fixture.dense, &model, &fixture.graph, &fixture.branch_lengths)?;
    let edge_key = fixture
      .graph
      .get_edges()
      .next()
      .map(|edge| edge.key())
      .expect("tree has edges");
    let plain = fixture
      .dense
      .marginal_update(&jc69, &fixture.graph, &fixture.branch_lengths, &())?;

    let single =
      fixture
        .dense
        .create_edge_contribution(&jc69, &plain.edges.backward, &plain.edges.forward, edge_key)?;
    let mixture = dense_mixture_contribution(&fixture.dense, &state, model.gamma().weights(), edge_key)?;

    assert_eq!(
      (true, None),
      (single.certified_decay_rate().is_some(), mixture.certified_decay_rate())
    );
    Ok(())
  }

  mod helpers {
    use crate::partition::optimize::contribution::OptimizationContribution;
    use crate::partition::storage::dense::{DenseEdgeEstimate, DenseNodeState};
    use crate::partition::storage::sparse::SparseNodeState;
    use crate::seq::indel::InDel;
    use crate::seq::mutation::Sub;
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{Array2, Axis, aview1, concatenate};
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::key_map::{EdgeMap, NodeMap};
    use treetime_primitives::Seq;

    pub(super) fn dense_profiles(nodes: &NodeMap<DenseNodeState>) -> Result<Array2<f64>, Report> {
      let views = nodes.values().map(|node| node.profile.dis.view()).collect_vec();
      Ok(concatenate(Axis(0), &views)?)
    }

    pub(super) fn dense_sequences(nodes: &NodeMap<DenseNodeState>) -> Vec<Seq> {
      nodes.values().map(|node| node.seq.sequence.clone()).collect()
    }

    pub(super) fn edge_indels(estimates: &EdgeMap<DenseEdgeEstimate>, edges: &[GraphEdgeKey]) -> Vec<Vec<InDel>> {
      edges
        .iter()
        .map(|edge_key| estimates[edge_key].indels.clone())
        .collect()
    }

    pub(super) fn sparse_estimates(estimates: &EdgeMap<Vec<Sub>>, edges: &[GraphEdgeKey]) -> Vec<Vec<Sub>> {
      edges.iter().map(|edge_key| estimates[edge_key].clone()).collect()
    }

    pub(super) fn stacked(contributions: &[OptimizationContribution]) -> Result<Array2<f64>, Report> {
      let views = contributions
        .iter()
        .filter_map(|contribution| match contribution {
          OptimizationContribution::Dense(contribution) => Some(contribution.coefficients.view()),
          OptimizationContribution::Sparse(_) => None,
        })
        .collect_vec();
      Ok(concatenate(Axis(0), &views)?)
    }

    pub(super) fn sparse_sites(contribution: &OptimizationContribution) -> Vec<(u64, Vec<u64>)> {
      match contribution {
        OptimizationContribution::Sparse(contribution) => contribution
          .site_contributions
          .iter()
          .map(|site| {
            (
              site.multiplicity.to_bits(),
              site.coefficients.iter().map(|value| value.to_bits()).collect_vec(),
            )
          })
          .sorted()
          .collect(),
        OptimizationContribution::Dense(_) => vec![],
      }
    }

    pub(super) fn sparse_posteriors(nodes: &NodeMap<SparseNodeState>, length: usize) -> Result<Array2<f64>, Report> {
      let rows = nodes
        .values()
        .flat_map(|node| (0..length).map(move |pos| posterior(node, pos)))
        .collect_vec();
      let views = rows.iter().map(|row| row.view().insert_axis(Axis(0))).collect_vec();
      Ok(concatenate(Axis(0), &views)?)
    }

    fn posterior(node: &SparseNodeState, pos: usize) -> ndarray::Array1<f64> {
      node.profile.variable.get(pos).map_or_else(
        || {
          node
            .profile
            .fixed
            .get(node.sequence[pos])
            .map_or_else(|| ndarray::Array1::zeros(4), |fixed| aview1(&fixed.dis).to_owned())
        },
        |var| aview1(var.dis).to_owned(),
      )
    }
  }
}
