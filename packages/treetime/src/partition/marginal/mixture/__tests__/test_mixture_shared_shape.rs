#[cfg(test)]
mod tests {
  use crate::gtr::gamma::{GammaRates, INITIAL_GAMMA_SHAPE};
  use crate::gtr::partition_model::PartitionModel;
  use crate::partition::marginal::mixture::__tests__::fixtures::tests::helpers::{
    GAP_FREE, LONG_BRANCHES, MISSING_DATA, scaled_fixture,
  };
  use crate::partition::marginal::reconstruction::DenseReconstruction;
  use crate::partition::set::PartitionSet;
  use crate::pretty_assert_ulps_eq;
  use crate::progress::NoopProgress;
  use crate::test_utils::make_nonuniform_gtr;
  use eyre::Report;
  use itertools::Itertools;
  use ndarray::array;

  #[test]
  fn test_partitions_with_identical_data_that_share_the_shape_get_the_shape_of_one_partition() -> Result<(), Report> {
    let fixture = scaled_fixture(&GAP_FREE, LONG_BRANCHES)?;
    let model = PartitionModel::new(make_nonuniform_gtr()?, GammaRates::new(4, INITIAL_GAMMA_SHAPE)?);
    let reconstruction =
      |shared: bool| DenseReconstruction::seeded(fixture.dense.clone(), model.clone().with_shared_shape(shared));

    let single = PartitionSet::Dense(vec![reconstruction(false)]).estimate_gamma_shapes(
      &fixture.graph,
      &fixture.branch_lengths,
      &NoopProgress,
    )?;
    let pair = PartitionSet::Dense(vec![reconstruction(true), reconstruction(true)]).estimate_gamma_shapes(
      &fixture.graph,
      &fixture.branch_lengths,
      &NoopProgress,
    )?;

    let shape = |set: &PartitionSet, index: usize| {
      set
        .iter()
        .nth(index)
        .and_then(|partition| partition.model().gamma().shape())
    };
    let expected = shape(&single, 0).expect("the single partition has a gamma shape");
    let actual = array![
      shape(&pair, 0).expect("the first partition has a gamma shape"),
      shape(&pair, 1).expect("the second partition has a gamma shape")
    ];
    pretty_assert_ulps_eq!(&array![expected, expected], &actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_partitions_that_share_the_shape_get_one_shape_and_others_their_own() -> Result<(), Report> {
    let first = scaled_fixture(&GAP_FREE, LONG_BRANCHES)?;
    let second = scaled_fixture(&MISSING_DATA, LONG_BRANCHES)?;
    let model = PartitionModel::new(make_nonuniform_gtr()?, GammaRates::new(4, INITIAL_GAMMA_SHAPE)?);
    let estimate = |shared: bool| -> Result<(f64, f64), Report> {
      let set = PartitionSet::Dense(vec![
        DenseReconstruction::seeded(first.dense.clone(), model.clone().with_shared_shape(shared)),
        DenseReconstruction::seeded(second.dense.clone(), model.clone().with_shared_shape(shared)),
      ])
      .estimate_gamma_shapes(&first.graph, &first.branch_lengths, &NoopProgress)?;
      let shapes = set
        .iter()
        .map(|partition| {
          partition
            .model()
            .gamma()
            .shape()
            .expect("every partition has a gamma shape")
        })
        .collect_vec();
      Ok((shapes[0], shapes[1]))
    };

    let (shared_first, shared_second) = estimate(true)?;
    let (own_first, own_second) = estimate(false)?;

    pretty_assert_ulps_eq!(shared_first, shared_second, epsilon = 1e-15);
    assert!(
      (own_first - own_second).abs() > 1e-3,
      "partitions with different data and their own shapes got {own_first} and {own_second}"
    );
    Ok(())
  }
}
