#[cfg(test)]
mod tests {
  use crate::array::ndarray::{argmax_all_within, argmax_first, argmax_first_within};
  use itertools::Itertools;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::single(               (array![0.3],                     1e-10), Some(0))]
  #[case::clear_maximum(        (array![0.2, 0.7, 0.1],           1e-10), Some(1))]
  #[case::exact_tie(            (array![0.1, 0.45, 0.45],         1e-10), Some(1))]
  #[case::later_larger_in_tie(  (array![0.5 - 4e-11, 0.5 + 4e-11], 1e-10), Some(0))]
  #[case::outside_tolerance(    (array![0.5 - 1e-9, 0.5 + 1e-9],  1e-10), Some(1))]
  #[case::zero_tolerance(       (array![0.4, 0.4 + 1e-15],        0.0),   Some(1))]
  #[case::nan_ignored(          (array![f64::NAN, 0.2, 0.8],      1e-10), Some(2))]
  #[case::all_nan(              (array![f64::NAN, f64::NAN],      1e-10), None)]
  #[case::empty(                (Array1::<f64>::zeros(0),         1e-10), None)]
  #[trace]
  fn test_argmax_first_within_takes_first_value_within_tolerance_of_maximum(
    #[case] (values, tolerance): (Array1<f64>, f64),
    #[case] expected: Option<usize>,
  ) {
    assert_eq!(expected, argmax_first_within(&values, tolerance));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_state(          (array![0.1, 0.8, 0.1],        0.125), vec![1])]
  #[case::two_states(         (array![0.45, 0.1, 0.45],      0.0),   vec![0, 2])]
  #[case::within_threshold(   (array![0.5, 0.4, 0.1],        0.125), vec![0, 1])]
  #[case::boundary_included(  (array![0.5, 0.25, 0.25],      0.25),  vec![0, 1, 2])]
  #[case::uniform(            (array![0.25, 0.25, 0.25, 0.25], 0.0), vec![0, 1, 2, 3])]
  #[case::empty(              (Array1::<f64>::zeros(0),      0.5),   vec![])]
  #[trace]
  fn test_argmax_all_within_lists_every_value_within_tolerance_of_maximum(
    #[case] (values, tolerance): (Array1<f64>, f64),
    #[case] expected: Vec<usize>,
  ) {
    assert_eq!(expected, argmax_all_within(&values, tolerance).collect_vec());
  }

  proptest! {
    #[test]
    fn test_prop_argmax_first_within_is_first_of_argmax_all_within(
      values in prop::collection::vec(0.0..1.0_f64, 1..12),
      tolerance in 0.0..0.5_f64,
    ) {
      let values = Array1::from_iter(values);
      prop_assert_eq!(argmax_all_within(&values, tolerance).next(), argmax_first_within(&values, tolerance));
    }

    #[test]
    fn test_prop_argmax_all_within_contains_the_maximum_and_grows_with_tolerance(
      values in prop::collection::vec(0.0..1.0_f64, 1..12),
      tolerance in 0.0..0.5_f64,
      extra in 0.0..0.5_f64,
    ) {
      let values = Array1::from_iter(values);
      let maximum = argmax_first(&values.view()).unwrap();
      let narrow = argmax_all_within(&values, tolerance).collect_vec();
      let wide = argmax_all_within(&values, tolerance + extra).collect_vec();
      prop_assert!(narrow.contains(&maximum));
      prop_assert!(narrow.iter().all(|index| wide.contains(index)));
    }
  }
}
