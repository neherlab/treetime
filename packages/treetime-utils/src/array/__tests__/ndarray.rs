#[cfg(test)]
mod tests {
  #![allow(clippy::excessive_precision, clippy::lossy_float_literal)]

  use std::sync::LazyLock;

  use crate::array::ndarray::*;
  use crate::pretty_assert_ulps_eq;
  use ::ndarray::{Array0, Array1, Array2, Axis, arr0, array};
  use eyre::Report;
  use rand::SeedableRng;
  use rand_isaac::Isaac64Rng;
  use rstest::rstest;

  static INPUT: LazyLock<Array2<f64>> = LazyLock::new(|| {
    array![
      [0.19356424, 0.25224431, 0.21259213, 0.19217803, 0.14942128],
      [0.19440831, 0.13170981, 0.26841564, 0.29005381, 0.11541244],
      [0.27439982, 0.18330691, 0.19687558, 0.32079767, 0.02462001],
      [0.03366488, 0.00781195, 0.32170632, 0.30066296, 0.33615390],
      [0.31185458, 0.25466645, 0.14705881, 0.24872985, 0.03769030],
      [0.24016971, 0.05380214, 0.35454510, 0.19585567, 0.15562739],
      [0.12705805, 0.37184099, 0.21907519, 0.27300161, 0.00902417],
    ]
  });

  #[rstest]
  fn computes_argmin_axis_0() {
    assert_eq!(argmin_axis(&INPUT, Axis(0)), array![3, 3, 4, 0, 6]);
  }

  #[rstest]
  fn computes_argmin_axis_1() {
    assert_eq!(argmin_axis(&INPUT, Axis(1)), array![4, 4, 4, 1, 4, 1, 4]);
  }

  #[rstest]
  fn computes_argmax_axis_0() {
    assert_eq!(argmax_axis(&INPUT, Axis(0)), array![4, 6, 5, 2, 3]);
  }

  #[rstest]
  fn computes_argmax_axis_1() {
    assert_eq!(argmax_axis(&INPUT, Axis(1)), array![1, 3, 3, 4, 0, 2, 1]);
  }

  #[rstest]
  fn computes_cumsum_axis_0() {
    let expected = array![
      [0.19356424, 0.25224431, 0.21259213, 0.19217803, 0.14942128],
      [0.38797255, 0.38395412, 0.48100777, 0.48223184, 0.26483372],
      [0.66237237, 0.56726103, 0.67788335, 0.80302951, 0.28945373],
      [0.69603725, 0.57507298, 0.99958967, 1.10369247, 0.62560763],
      [1.00789183, 0.82973943, 1.14664848, 1.35242232, 0.66329793],
      [1.24806154, 0.88354157, 1.50119358, 1.54827799, 0.81892532],
      [1.37511959, 1.25538256, 1.72026877, 1.82127960, 0.82794949],
    ];

    pretty_assert_ulps_eq!(cumsum_axis(&INPUT, Axis(0)), expected);
  }

  #[rstest]
  fn computes_cumsum_axis_1() {
    let expected = array![
      [0.19356424, 0.44580855, 0.65840068, 0.85057871, 0.99999999],
      [0.19440831, 0.32611812, 0.59453376, 0.88458757, 1.00000001],
      [0.27439982, 0.45770673, 0.65458231, 0.97537998, 0.99999999],
      [0.03366488, 0.04147683, 0.36318315, 0.66384611, 1.00000001],
      [0.31185458, 0.56652103, 0.71357984, 0.96230969, 0.99999999],
      [0.24016971, 0.29397185, 0.64851695, 0.84437262, 1.00000001],
      [0.12705805, 0.49889904, 0.71797423, 0.99097584, 1.00000001],
    ];

    pretty_assert_ulps_eq!(cumsum_axis(&INPUT, Axis(1)), expected);
  }

  #[rstest]
  fn computes_outer_product() -> Result<(), Report> {
    pretty_assert_ulps_eq!(
      outer(&array![0.0, 1.0, 2.0, 3.0, 4.0], &array![-2.0, -1.0, 0.0, 1.0, 2.0])?,
      array![
        [-0.0, -0.0, 0.0, 0.0, 0.0],
        [-2.0, -1.0, 0.0, 1.0, 2.0],
        [-4.0, -2.0, 0.0, 2.0, 4.0],
        [-6.0, -3.0, 0.0, 3.0, 6.0],
        [-8.0, -4.0, 0.0, 4.0, 8.0]
      ]
    );
    Ok(())
  }

  #[test]
  fn test_outer_product_case_1() {
    let pi = array![0.25, 0.25, 0.25, 0.25];
    let ti = array![1.98000002, 2.94500003, 2.51500002, 2.64000002];
    #[rustfmt::skip]
    let expected_outer_pi_ti = array![
      [0.4950000050, 0.7362500075, 0.6287500050, 0.6600000050],
      [0.4950000050, 0.7362500075, 0.6287500050, 0.6600000050],
      [0.4950000050, 0.7362500075, 0.6287500050, 0.6600000050],
      [0.4950000050, 0.7362500075, 0.6287500050, 0.6600000050],
    ];
    #[rustfmt::skip]
    let expected_outer_ti_pi = array![
      [0.4950000050, 0.4950000050, 0.4950000050, 0.4950000050],
      [0.7362500075, 0.7362500075, 0.7362500075, 0.7362500075],
      [0.6287500050, 0.6287500050, 0.6287500050, 0.6287500050],
      [0.6600000050, 0.6600000050, 0.6600000050, 0.6600000050],
    ];
    let result_outer_pi_ti = outer(&pi, &ti).unwrap();
    let result_outer_ti_pi = outer(&ti, &pi).unwrap();
    pretty_assert_ulps_eq!(expected_outer_pi_ti, result_outer_pi_ti, max_ulps = 4);
    pretty_assert_ulps_eq!(expected_outer_ti_pi, result_outer_ti_pi, max_ulps = 4);
  }

  #[test]
  fn test_outer_product_case_2() {
    let pi = array![0.14873251, 0.24050277, 0.31232323, 0.29844148];
    let ti = array![1.98000002, 2.94500003, 2.51500002, 2.64000002];
    #[rustfmt::skip]
    let expected_outer_pi_ti = array![
      [0.29449037277465023, 0.43801724641197537, 0.37406226562465020, 0.39265382937465020],
      [0.47619548941005540, 0.70828066486508310, 0.60486447136005540, 0.63492731761005540],
      [0.61840000164646470, 0.91979192171969690, 0.78549292969646470, 0.82453333344646460],
      [0.59091413636882960, 0.87891016755324440, 0.75058032816882960, 0.78788551316882960],
    ];
    #[rustfmt::skip]
    let expected_outer_ti_pi = array![
      [0.29449037277465023, 0.47619548941005540, 0.61840000164646470, 0.59091413636882960],
      [0.43801724641197537, 0.70828066486508310, 0.91979192171969690, 0.87891016755324440],
      [0.37406226562465020, 0.60486447136005540, 0.78549292969646470, 0.75058032816882960],
      [0.39265382937465020, 0.63492731761005540, 0.82453333344646460, 0.78788551316882960],
    ];
    let result_outer_pi_ti = outer(&pi, &ti).unwrap();
    let result_outer_ti_pi = outer(&ti, &pi).unwrap();
    pretty_assert_ulps_eq!(expected_outer_pi_ti, result_outer_pi_ti, max_ulps = 4);
    pretty_assert_ulps_eq!(expected_outer_ti_pi, result_outer_ti_pi, max_ulps = 4);
  }

  #[rstest]
  fn chooses_1d_by_index() {
    let indices = array![2, 0, 1, 0, 0, 2];
    let choices = array!['a', 'b', 'c'];
    assert_eq!(choose1(&indices, &choices), array!['c', 'a', 'b', 'a', 'a', 'c']);
  }

  #[rstest]
  fn chooses_2d_by_index() {
    let indices = array![2, 0, 1];
    let choices = array![[0, 1, 2], [10, 11, 12], [20, 21, 22], [30, 31, 32]];
    assert_eq!(choose2(&indices, &choices), array![20, 1, 12]);
  }

  #[rstest]
  fn generates_predictable_random_uniform() {
    let mut rng = Isaac64Rng::seed_from_u64(42);

    let r: Array2<f64> = random((3, 4), &mut rng);
    #[rustfmt::skip]
    let expected = array![
      [0.47098793460940813, 0.77817923328818700, 0.72803987044229210, 0.24949100191426288],
      [0.41172712062120720, 0.76245969702116350, 0.79208560304293950, 0.68138502792473000],
      [0.67137133889579470, 0.85151789643141470, 0.31181959423520780, 0.73982541135527980],
    ];
    pretty_assert_ulps_eq!(r, expected);

    let r: Array2<f64> = random((3, 4), &mut rng);
    #[rustfmt::skip]
    let expected = array![
      [0.09054507828914282, 0.72883120582400560, 0.04697932861607845, 0.58466146369624310],
      [0.60534170114191200, 0.23101797776198140, 0.66710289504536500, 0.92825277017181910],
      [0.07155806372486850, 0.21012760390829310, 0.60105776255346570, 0.94540503820318000],
    ];
    pretty_assert_ulps_eq!(r, expected);
  }

  #[test]
  fn test_product_axis_empty_array1() {
    let input: Array1<f64> = array![];
    let expected: Array0<f64> = arr0(1.0);
    let result = product_axis(&input, Axis(0));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_product_axis_empty_array2_axis0() {
    let input: Array2<f64> = array![[]];
    let expected: Array1<f64> = array![];
    let result = product_axis(&input, Axis(0));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_product_axis_empty_array2_axis1() {
    let input: Array2<f64> = array![[]];
    let expected: Array1<f64> = array![1.0];
    let result = product_axis(&input, Axis(1));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_product_axis_empty_axis_0_prod_axis_0() {
    let input: Array2<f64> = array![[2.0, 3.0, 4.0]];
    let expected: Array1<f64> = array![2.0, 3.0, 4.0];
    let result = product_axis(&input, Axis(0));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_product_axis_empty_axis_0_prod_axis_1() {
    let input: Array2<f64> = array![[2.0, 3.0, 4.0]];
    let expected: Array1<f64> = array![24.0];
    let result = product_axis(&input, Axis(1));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_product_axis_empty_axis_1_prod_axis_0() {
    let input: Array2<f64> = array![[2.0], [3.0], [4.0]];
    let expected: Array1<f64> = array![24.0];
    let result = product_axis(&input, Axis(0));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_product_axis_empty_axis_1_prod_axis_1() {
    let input: Array2<f64> = array![[2.0], [3.0], [4.0]];
    let expected: Array1<f64> = array![2.0, 3.0, 4.0];
    let result = product_axis(&input, Axis(1));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_product_axis_general_case_axis_0() {
    let input = array![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]];
    let expected = array![15.0, 48.0];
    let result = product_axis(&input, Axis(0));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }

  #[test]
  fn test_product_axis_general_case_axis_1() {
    let input = array![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]];
    let expected = array![2.0, 12.0, 30.0];
    let result = product_axis(&input, Axis(1));
    pretty_assert_ulps_eq!(result, expected, max_ulps = 4);
  }
}
