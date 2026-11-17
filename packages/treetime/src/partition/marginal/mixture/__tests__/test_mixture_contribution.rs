#[cfg(test)]
mod tests {
  use crate::optimize::dispatch::optimize_edge_lengths;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::likelihood::evaluate_partitions_log_lh_only;
  use crate::optimize::params::BranchOptMethod;
  use crate::partition::marginal::mixture::__tests__::fixtures::tests::helpers::{
    GAP_FREE, LONG_BRANCHES, MISSING_DATA,
  };
  use crate::pretty_assert_abs_diff_eq;
  use eyre::Report;
  use helpers::{Representation, mixture_set};
  use itertools::Itertools;
  use ndarray::Array1;
  use rstest::rstest;
  use treetime_graph::key_map::EdgeMap;

  #[rustfmt::skip]
  #[rstest]
  #[case::dense_gap_free(  Representation::Dense,  &GAP_FREE)]
  #[case::sparse_gap_free( Representation::Sparse, &GAP_FREE)]
  #[case::dense_missing(   Representation::Dense,  &MISSING_DATA)]
  #[case::sparse_missing(  Representation::Sparse, &MISSING_DATA)]
  #[trace]
  fn test_mixture_contribution_differences_match_full_passes(
    #[case] representation: Representation,
    #[case] alignment: &[(&str, &str)],
  ) -> Result<(), Report> {
    let (fixture, set) = mixture_set(representation, alignment, LONG_BRANCHES, 5.0)?;
    let partition = set.iter().next().expect("set has one partition");
    let full = |edge_key, length: f64| -> Result<f64, Report> {
      let mut lengths = fixture.branch_lengths.clone();
      lengths.insert(edge_key, length);
      Ok(set.clone().marginal_update(&fixture.graph, &lengths)?.1.value())
    };
    let edges = fixture.graph.get_edges().map(|edge| edge.key()).collect_vec();

    let expected = edges
      .iter()
      .map(|&edge_key| Ok(full(edge_key, 0.4)? - full(edge_key, 0.05)?))
      .collect::<Result<Array1<f64>, Report>>()?;
    let actual = edges
      .iter()
      .map(|&edge_key| {
        let contribution = [partition.create_edge_contribution(edge_key)?];
        Ok(
          evaluate_partitions_log_lh_only(&contribution, 0.4)?.value()
            - evaluate_partitions_log_lh_only(&contribution, 0.05)?.value(),
        )
      })
      .collect::<Result<Array1<f64>, Report>>()?;

    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense_brent(  Representation::Dense,  BranchOptMethod::Brent)]
  #[case::dense_newton( Representation::Dense,  BranchOptMethod::Newton)]
  #[case::sparse_brent( Representation::Sparse, BranchOptMethod::Brent)]
  #[case::sparse_newton(Representation::Sparse, BranchOptMethod::Newton)]
  #[trace]
  fn test_branch_length_methods_find_the_optimum_of_the_full_mixture_objective(
    #[case] representation: Representation,
    #[case] method: BranchOptMethod,
  ) -> Result<(), Report> {
    let (fixture, set) = mixture_set(representation, &GAP_FREE, 1.0, 0.5)?;
    let partition = set.iter().next().expect("set has one partition");
    let mut lengths: EdgeMap<Option<f64>> = fixture.branch_lengths.iter().map(|(key, &length)| (key, Some(length))).collect();
    let indels = IndelPolicy::Exclude.edge_indels(&fixture.graph, &set, &lengths)?;
    optimize_edge_lengths(&fixture.graph, set.total_length(), &set, &indels, method, &mut lengths)?;

    let root = fixture.graph.root_key()?;
    let full = |edge_key, length: f64| -> Result<f64, Report> {
      let mut other = fixture.branch_lengths.clone();
      other.insert(edge_key, length);
      Ok(set.clone().marginal_update(&fixture.graph, &other)?.1.value())
    };
    let gains = fixture
      .graph
      .get_edges()
      .map(|edge| edge.key())
      .filter(|&edge_key| fixture.graph.edge_endpoints(edge_key).is_ok_and(|(parent, _)| parent != root))
      .map(|edge_key| {
        let optimum = lengths[&edge_key].expect("optimized length");
        let neighbour = full(edge_key, optimum * 0.97)?.max(full(edge_key, optimum * 1.03 + 1e-6)?);
        Ok(full(edge_key, optimum)? - neighbour)
      })
      .collect::<Result<Vec<f64>, Report>>()?;

    let worst = gains.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(worst >= -1e-9, "an optimized length is no optimum of the full mixture objective: {gains:?}");
    Ok(())
  }

  mod helpers {
    pub(super) const CATEGORIES: usize = 4;
    use crate::gtr::gamma::GammaRates;
    use crate::gtr::partition_model::PartitionModel;
    use crate::partition::marginal::mixture::__tests__::fixtures::tests::helpers::{Fixture, scaled_fixture};
    use crate::partition::marginal::reconstruction::{DenseReconstruction, SparseReconstruction};
    use crate::partition::set::PartitionSet;
    use crate::test_utils::make_nonuniform_gtr;
    use eyre::Report;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Representation {
      Dense,
      Sparse,
    }

    pub(super) fn mixture_set(
      representation: Representation,
      alignment: &[(&str, &str)],
      scale: f64,
      shape: f64,
    ) -> Result<(Fixture, PartitionSet), Report> {
      let fixture = scaled_fixture(alignment, scale)?;
      let model = PartitionModel::new(make_nonuniform_gtr()?, GammaRates::new(CATEGORIES, shape)?);
      let set = match representation {
        Representation::Dense => PartitionSet::Dense(vec![
          DenseReconstruction::seeded(fixture.dense.clone(), model)
            .marginal_update(&fixture.graph, &fixture.branch_lengths)?
            .0,
        ]),
        Representation::Sparse => PartitionSet::Sparse(vec![
          SparseReconstruction::seeded(fixture.sparse.clone(), model, fixture.sparse_nodes.clone())
            .marginal_update(&fixture.graph, &fixture.branch_lengths)?
            .0,
        ]),
      };
      Ok((fixture, set))
    }
  }
}
