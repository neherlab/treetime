#[cfg(test)]
mod tests {
  use crate::gtr::gamma::GammaRates;
  use eyre::Report;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

  #[test]
  fn test_gamma_rates_with_one_category_are_the_single_rate_one() -> Result<(), Report> {
    let expected = GammaRates::default();
    let actual = GammaRates::new(1, 0.5)?;
    assert_eq!(expected, actual);
    assert_eq!(&array![1.0], actual.rates());
    assert_eq!(None, actual.shape());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_categories(0, 1.0,   "Number of gamma rate categories must be at least 1, got 0")]
  #[case::shape_below(  4, 0.01,  "Gamma shape must lie between 0.02 and 100, got 0.01")]
  #[case::shape_above(  4, 101.0, "Gamma shape must lie between 0.02 and 100, got 101")]
  #[case::shape_nan(    4, f64::NAN, "Gamma shape must lie between 0.02 and 100, got NaN")]
  #[trace]
  fn test_gamma_rates_reject_invalid_arguments(#[case] categories: usize, #[case] shape: f64, #[case] expected: &str) {
    assert_error!(GammaRates::new(categories, shape), expected);
  }

  #[test]
  fn test_gamma_rates_increase_with_the_category() -> Result<(), Report> {
    let rates = GammaRates::new(4, 0.5)?;
    assert!(
      rates.rates().windows(2).into_iter().all(|pair| pair[0] < pair[1]),
      "the mean of a higher quantile interval is larger, got {}",
      rates.rates()
    );
    assert_eq!(&array![0.25, 0.25, 0.25, 0.25], rates.weights());
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_gamma_rates_have_mean_one(categories in 2_usize..12, shape in 0.02_f64..100.0) {
      let rates = GammaRates::new(categories, shape).unwrap();
      let mean = rates.rates().dot(rates.weights());
      pretty_assert_ulps_eq!(1.0, mean, epsilon = 1e-12);
      prop_assert!(rates.rates().iter().all(|rate| rate.is_finite() && *rate >= 0.0));
    }
  }
}
