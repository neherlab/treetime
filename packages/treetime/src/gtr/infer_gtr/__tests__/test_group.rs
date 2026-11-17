#[cfg(test)]
mod tests {
  use crate::gtr::gtr::{GTR, avg_transition};
  use crate::gtr::infer_gtr::common::{
    GtrSharing, InferGtrOptions, InferGtrResult, MutationCounts, infer_gtr_from_counts, infer_gtr_group,
  };
  use crate::pretty_assert_ulps_eq;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use helpers::{changing, changing_other, offdiagonal_ratios, stored_gtr, without_changes};
  use ndarray::{Array1, array};
  use pretty_assertions::assert_ne;

  const SHARE_ALL: GtrSharing = GtrSharing {
    exchange_rates: true,
    frequencies: true,
  };

  #[test]
  fn test_identical_partitions_that_share_everything_fit_like_their_pooled_counts() -> Result<(), Report> {
    let options = InferGtrOptions {
      pc: 0.0,
      ..InferGtrOptions::default()
    };
    let counts = changing();
    let pooled = MutationCounts {
      nij: 2.0 * &counts.nij,
      Ti: 2.0 * &counts.Ti,
      root_state: 2.0 * &counts.root_state,
    };

    let group = infer_gtr_group(&[&counts, &counts], SHARE_ALL, &options, &NoopProgress)?;
    let expected = infer_gtr_from_counts(&pooled, &options, &NoopProgress)?;

    pretty_assert_ulps_eq!(&expected.W, &group[0].W, epsilon = 1e-12);
    pretty_assert_ulps_eq!(&expected.pi, &group[0].pi, epsilon = 1e-12);
    pretty_assert_ulps_eq!(
      &array![expected.mu.value(), expected.mu.value()],
      &array![group[0].mu.value(), group[1].mu.value()],
      epsilon = 1e-12
    );
    pretty_assert_ulps_eq!(&group[0].W, &group[1].W, epsilon = 1e-15);
    pretty_assert_ulps_eq!(&group[0].pi, &group[1].pi, epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_shared_fit_adds_the_pseudocount_once() -> Result<(), Report> {
    let counts = changing();
    let empty = MutationCounts {
      nij: 0.0 * &counts.nij,
      Ti: 0.0 * &counts.Ti,
      root_state: 0.0 * &counts.root_state,
    };

    let group = infer_gtr_group(
      &[&counts, &empty],
      SHARE_ALL,
      &InferGtrOptions::default(),
      &NoopProgress,
    )?;
    let expected = infer_gtr_from_counts(&counts, &InferGtrOptions::default(), &NoopProgress)?;

    pretty_assert_ulps_eq!(&expected.W, &group[0].W, epsilon = 1e-12);
    pretty_assert_ulps_eq!(&expected.pi, &group[0].pi, epsilon = 1e-12);
    pretty_assert_ulps_eq!(expected.mu.value(), group[0].mu.value(), epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_shared_exchange_rates_differ_only_by_a_scale_under_different_frequencies() -> Result<(), Report> {
    let sharing = GtrSharing {
      exchange_rates: true,
      frequencies: false,
    };
    let (first, second) = (changing(), changing_other());

    let group = infer_gtr_group(&[&first, &second], sharing, &InferGtrOptions::default(), &NoopProgress)?;
    let gtrs = group.iter().map(stored_gtr).collect::<Result<Vec<_>, Report>>()?;

    assert_ne!(gtrs[0].pi(), gtrs[1].pi());
    let ratios = offdiagonal_ratios(gtrs[0].W(), gtrs[1].W());
    let scale = ratios[0];
    pretty_assert_ulps_eq!(&Array1::from_elem(ratios.len(), scale), &ratios, epsilon = 1e-12);
    assert!(
      (scale - 1.0).abs() > 1e-3,
      "the stored exchange rates differ by the factor {scale}"
    );
    let mean_rates = gtrs
      .iter()
      .map(|gtr| avg_transition(gtr.W(), gtr.pi()))
      .collect::<Array1<f64>>();
    pretty_assert_ulps_eq!(&array![1.0, 1.0], &mean_rates, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_partition_without_changes_takes_the_shared_exchange_rates_and_keeps_its_rate() -> Result<(), Report> {
    let sharing = GtrSharing {
      exchange_rates: true,
      frequencies: false,
    };
    let (first, quiet) = (changing(), without_changes());

    let group = infer_gtr_group(&[&first, &quiet], sharing, &InferGtrOptions::default(), &NoopProgress)?;
    let gtrs = group.iter().map(stored_gtr).collect::<Result<Vec<_>, Report>>()?;

    let ratios = offdiagonal_ratios(gtrs[0].W(), gtrs[1].W());
    pretty_assert_ulps_eq!(&Array1::from_elem(ratios.len(), ratios[0]), &ratios, epsilon = 1e-12);
    assert!(
      gtrs[1].mu() < gtrs[0].mu(),
      "the partition without changes has rate {} and the changing partition {}",
      gtrs[1].mu(),
      gtrs[0].mu()
    );
    Ok(())
  }

  #[test]
  fn test_shared_frequencies_keep_one_set_of_exchange_rates_per_partition() -> Result<(), Report> {
    let sharing = GtrSharing {
      exchange_rates: false,
      frequencies: true,
    };
    let (first, second) = (changing(), changing_other());

    let group = infer_gtr_group(&[&first, &second], sharing, &InferGtrOptions::default(), &NoopProgress)?;

    pretty_assert_ulps_eq!(&group[0].pi, &group[1].pi, epsilon = 1e-15);
    let ratios = offdiagonal_ratios(&group[0].W, &group[1].W);
    assert!(
      (ratios[0] - ratios[1]).abs() > 1e-3,
      "exchange rates of partitions with their own fits are not proportional: {ratios}"
    );
    Ok(())
  }

  mod helpers {
    use super::{GTR, InferGtrResult, MutationCounts};
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{Array1, Array2, array};

    pub(super) fn changing() -> MutationCounts {
      MutationCounts {
        nij: array![
          [0.0, 1.0, 2.0, 1.0],
          [1.0, 0.0, 3.0, 2.0],
          [2.0, 3.0, 0.0, 1.0],
          [2.0, 3.0, 3.0, 0.0]
        ],
        Ti: array![12.0, 20.0, 14.0, 12.4],
        root_state: array![3.0, 2.0, 3.0, 4.0],
      }
    }

    pub(super) fn changing_other() -> MutationCounts {
      MutationCounts {
        nij: array![
          [0.0, 6.0, 1.0, 0.0],
          [4.0, 0.0, 1.0, 5.0],
          [1.0, 0.0, 0.0, 7.0],
          [0.0, 4.0, 6.0, 0.0]
        ],
        Ti: array![30.0, 8.0, 9.0, 21.0],
        root_state: array![6.0, 1.0, 1.0, 4.0],
      }
    }

    pub(super) fn without_changes() -> MutationCounts {
      MutationCounts {
        nij: Array2::zeros((4, 4)),
        Ti: array![10.0, 10.0, 10.0, 10.0],
        root_state: array![2.0, 2.0, 2.0, 2.0],
      }
    }

    pub(super) fn stored_gtr(result: &InferGtrResult) -> Result<GTR, Report> {
      GTR::builder()
        .n_states(4)
        .mu(result.mu)
        .W(result.W.clone())
        .pi(result.pi.clone())
        .build()
    }

    pub(super) fn offdiagonal_ratios(first: &Array2<f64>, second: &Array2<f64>) -> Array1<f64> {
      (0..4)
        .cartesian_product(0..4)
        .filter(|(i, j)| i != j)
        .map(|(i, j)| first[[i, j]] / second[[i, j]])
        .collect()
    }
  }
}
