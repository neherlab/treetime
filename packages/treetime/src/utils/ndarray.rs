use crate::constants::BIG_NUMBER;
use eyre::Report;
use itertools::Itertools;
use ndarray::{
  stack, Array, Array1, Array2, ArrayBase, ArrayView, Axis, Data, Dimension, Ix1, Ix2, RemoveAxis, ScalarOperand,
  ShapeBuilder, ShapeError, Zip,
};
use ndarray_rand::RandomExt;
use num_traits::real::Real;
use num_traits::{Bounded, Float, NumCast, One, Zero};
use rand::distributions::uniform::SampleUniform;
use rand::distributions::Uniform;
use rand::Rng;
use std::f64::consts::E;
use std::ops::AddAssign;

pub fn to_col<T: Real>(a: &Array1<T>) -> Result<Array2<T>, Report> {
  Ok(a.to_shape((a.len(), 1))?.into_dimensionality::<Ix2>()?.to_owned())
}

pub fn to_row<T: Real>(b: &Array1<T>) -> Result<Array2<T>, Report> {
  Ok(b.to_shape((1, b.len()))?.into_dimensionality::<Ix2>()?.to_owned())
}

// Stack owned arrays. Similar to ndarray::stack() but accepting owned arrays rather than views.
pub fn stack_owned<A, D>(axis: Axis, arrays: &[Array<A, D>]) -> Result<Array<A, D::Larger>, ShapeError>
where
  A: Clone + Default,
  D: Dimension,
  D::Larger: RemoveAxis,
{
  if arrays.is_empty() {
    return Ok(Array::<A, D::Larger>::default(D::Larger::default()));
  }
  let arrays = arrays.iter().map(ArrayView::from).collect_vec();
  stack(axis, &arrays)
}

// Calculates outer product of 2 vectors
pub fn outer<T: 'static + Real>(a: &Array1<T>, b: &Array1<T>) -> Result<Array2<T>, Report> {
  let a = a.to_shape((a.len(), 1))?.into_dimensionality::<Ix2>()?;
  let b = b.to_shape((1, b.len()))?.into_dimensionality::<Ix2>()?;
  Ok(a.dot(&b))
}

/// Calculates min over given axis
#[inline]
pub fn min_axis<T: Real>(arr: &Array2<T>, axis: Axis) -> Array1<T> {
  arr.fold_axis(axis, T::max_value(), |&a, &b| a.min(b))
}

/// Calculates max over given axis
#[inline]
pub fn max_axis<S, T: Real>(arr: &ArrayBase<S, Ix2>, axis: Axis) -> Array1<T>
where
  S: Data<Elem = T>,
{
  arr.fold_axis(axis, T::min_value(), |&a, &b| a.max(b))
}

/// Finds index of min value over given axis
#[inline]
pub fn argmin_axis<T: 'static + Real, D: RemoveAxis>(arr: &Array<T, D>, axis: Axis) -> Array<usize, D::Smaller> {
  arr
    .fold_axis(axis, (0_usize, 0_usize, T::max_value()), |(i_curr, i_min, x_min), x| {
      if x < x_min {
        (i_curr + 1, *i_curr, *x)
      } else {
        (i_curr + 1, *i_min, *x_min)
      }
    })
    .mapv_into_any(|(_, i, _)| i)
}

/// Finds index of max value over given axis
#[inline]
pub fn argmax_axis<T: 'static + Copy + PartialOrd + Bounded, D: RemoveAxis>(
  arr: &Array<T, D>,
  axis: Axis,
) -> Array<usize, D::Smaller> {
  arr
    .fold_axis(axis, (0_usize, 0_usize, T::min_value()), |(i_curr, i_max, x_max), x| {
      if x > x_max {
        (i_curr + 1, *i_curr, *x)
      } else {
        (i_curr + 1, *i_max, *x_max)
      }
    })
    .mapv_into_any(|(_, i, _)| i)
}

/// Element-wise minimum of two arrays
pub fn minimum<T: Copy + PartialOrd, D: Dimension>(x: &Array<T, D>, y: &Array<T, D>) -> Array<T, D> {
  assert_eq!(x.shape(), y.shape());
  Zip::from(x).and(y).map_collect(|&a, &b| if a < b { a } else { b })
}

/// Element-wise maximum of two arrays
pub fn maximum<T: Copy + PartialOrd, D: Dimension>(x: &Array<T, D>, y: &Array<T, D>) -> Array<T, D> {
  assert_eq!(x.shape(), y.shape());
  Zip::from(x).and(y).map_collect(|&a, &b| if a > b { a } else { b })
}

/// Element-wise maximum of array and a scalar
pub fn maximum_scalar<T: Copy + PartialOrd, D: Dimension>(arr: &Array<T, D>, y: T) -> Array<T, D> {
  arr.mapv(|a| if a > y { a } else { y })
}

/// Clamp each element to at most `lower`
pub fn clamp_min<T: Copy + PartialOrd, D: Dimension>(a: &Array<T, D>, lower: T) -> Array<T, D> {
  a.mapv(|x| num_traits::clamp_min(x, lower))
}

/// Clamp each element to at most `upper`
pub fn clamp_max<T: Copy + PartialOrd, D: Dimension>(a: &Array<T, D>, upper: T) -> Array<T, D> {
  a.mapv(|x| num_traits::clamp_max(x, upper))
}

/// Clamp each element so that they are between given `lower` and `upper` values
pub fn clamp<T: Copy + PartialOrd, D: Dimension>(a: &Array<T, D>, lower: T, upper: T) -> Array<T, D> {
  a.mapv(|x| num_traits::clamp(x, lower, upper))
}

/// Element-wise exp
pub fn exp<T, S, D>(arr: &ArrayBase<S, D>) -> Array<T, D>
where
  T: Float + From<f64>,
  S: Data<Elem = T>,
  D: Dimension,
{
  arr.mapv(|x| Float::exp(x))
}

/// Element-wise log
pub fn log<T: Float + From<f64>, D: Dimension>(arr: &Array<T, D>) -> Array<T, D> {
  arr.mapv(|x| Float::log(x, E.into()))
}

/// Construct an array from an index array and an array to choose from
pub fn choose1<T, SI, ST>(indices: &ArrayBase<SI, Ix1>, arr: &ArrayBase<ST, Ix1>) -> Array1<T>
where
  T: Copy + Default,
  SI: Data<Elem = usize>,
  ST: Data<Elem = T>,
{
  indices.iter().map(|i| arr[*i]).collect()
}

/// Construct an array from an index array and a list of arrays to choose from (or 2D array).
pub fn choose2<T, SI, ST>(indices: &ArrayBase<SI, Ix1>, arr: &ArrayBase<ST, Ix2>) -> Array1<T>
where
  T: Copy + Default,
  SI: Data<Elem = usize>,
  ST: Data<Elem = T>,
{
  Zip::from(indices)
    .and(arr.axis_iter(Axis(1)))
    .map_collect(|i, axis| axis[*i])
}

pub fn zeros<T, D, Sh>(shape: Sh) -> Array<T, D>
where
  T: Zero + Clone,
  D: Dimension,
  Sh: ShapeBuilder<Dim = D>,
{
  Array::<T, D>::zeros(shape)
}

/// Calculates cumulative sum over given axis
#[inline]
pub fn cumsum_axis<T: Copy + AddAssign, D: Dimension>(a: &Array<T, D>, axis: Axis) -> Array<T, D> {
  let mut result = a.to_owned();
  result.accumulate_axis_inplace(axis, |&prev, curr| *curr += prev);
  result
}

/// Calculate product over given axis
#[inline]
pub fn product_axis<T: Copy + One + ScalarOperand, D: Dimension + RemoveAxis>(
  a: &Array<T, D>,
  axis: Axis,
) -> Array<T, D::Smaller> {
  let init_value = Array::<T, D::Smaller>::ones(a.raw_dim().remove_axis(axis));
  a.fold_axis(axis, T::one(), |&acc, &x| acc * x)
}

pub fn random<T: Copy + SampleUniform + NumCast, D: Dimension, Sh: ShapeBuilder<Dim = D>, R: Rng>(
  shape: Sh,
  rng: &mut R,
) -> Array<T, D> {
  let from: T = NumCast::from(0_i32).unwrap();
  let to: T = NumCast::from(1_i32).unwrap();
  Array::<T, D>::random_using(shape, Uniform::<T>::new::<T, T>(from, to), rng)
}

/// Replaces non-finite with a very small number
pub fn sanitize_in_place<T: Float, D: Dimension>(arr: &mut Array<T, D>) {
  arr.mapv_inplace(|x| {
    if x.is_finite() {
      x
    } else {
      T::from(-BIG_NUMBER).unwrap()
    }
  });
}

#[cfg(test)]
mod tests {
  #![allow(clippy::excessive_precision, clippy::lossy_float_literal)]

  use super::*;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use lazy_static::lazy_static;
  use ndarray::array;
  use ndarray_linalg::{Eigh, UPLO};
  use rand::SeedableRng;
  use rand_isaac::Isaac64Rng;
  use rstest::rstest;

  lazy_static! {
    static ref INPUT: Array2<f64> = array![
      [0.19356424, 0.25224431, 0.21259213, 0.19217803, 0.14942128],
      [0.19440831, 0.13170981, 0.26841564, 0.29005381, 0.11541244],
      [0.27439982, 0.18330691, 0.19687558, 0.32079767, 0.02462001],
      [0.03366488, 0.00781195, 0.32170632, 0.30066296, 0.33615390],
      [0.31185458, 0.25466645, 0.14705881, 0.24872985, 0.03769030],
      [0.24016971, 0.05380214, 0.35454510, 0.19585567, 0.15562739],
      [0.12705805, 0.37184099, 0.21907519, 0.27300161, 0.00902417],
    ];
  }

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

  #[rstest]
  fn computes_eigh() -> Result<(), Report> {
    // Comparison of Rust ndarray_linalg::eigh() and NumPy np.linalg.eigh()
    // https://docs.rs/ndarray-linalg/latest/ndarray_linalg/eigh/index.html
    // https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html

    let a: Array2<f64> = array![
      [-1.0, 0.25, 0.25, 0.25, 0.25],
      [0.25, -1.0, 0.25, 0.25, 0.25],
      [0.25, 0.25, -1.0, 0.25, 0.25],
      [0.25, 0.25, 0.25, -1.0, 0.25],
      [0.25, 0.25, 0.25, 0.25, -1.0],
    ];

    // Rust version:
    let (eigvals, eigvecs) = a.eigh(UPLO::Lower)?;

    // NumPy version:
    // import numpy as np
    // np.set_printoptions(precision=60, suppress=True, linewidth=999, sign=' ', floatmode='maxprec_equal')
    // pprint(np.linalg.eigh([
    //   [-1.0, 0.25, 0.25, 0.25, 0.25],
    //   [0.25, -1.0, 0.25, 0.25, 0.25],
    //   [0.25, 0.25, -1.0, 0.25, 0.25],
    //   [0.25, 0.25, 0.25, -1.0, 0.25],
    //   [0.25, 0.25, 0.25, 0.25, -1.0],
    // ]))

    // WolframAlpha version:
    // https://www.wolframalpha.com/input?i2d=true&i=%7B%7B-1.0%2C0.25%2C0.25%2C0.25%2C0.25%7D%2C%7B0.25%2C-1.0%2C0.25%2C0.25%2C0.25%7D%2C%7B0.25%2C0.25%2C-1.0%2C0.25%2C0.25%7D%2C%7B0.25%2C0.25%2C0.25%2C-1.0%2C0.25%7D%2C%7B0.25%2C0.25%2C0.25%2C0.25%2C-1.0%7D%7D

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      eigvals,

      // NumPy np.linalg.eigh() result (with openblas backend):
      array![-1.25000000000000022204460492503131, -1.25000000000000000000000000000000, -1.24999999999999977795539507496869, -1.24999999999999888977697537484346,  0.00000000000000005551115123125783],

      max_ulps = 1
    );

    #[rustfmt::skip]
    pretty_assert_ulps_eq!(
      eigvecs,

      // NumPy np.linalg.eigh() result (with openblas backend):
      array![
        [ 0.00000000000000000000000000000000,  0.67791921943923838522749747426133,  0.58345996599157845530214672180591, -0.00000000000000000000000000000000,  0.44721359549995798321475604097941],
        [-0.00000000000000004163336342344337,  0.39545287800622275220518986316165, -0.80225745323841990419566627679160,  0.00000000000000004163336342344337,  0.44721359549995809423705850349506],
        [-0.28307239576681808568281439875136, -0.35779069914848726785550070417230,  0.07293249574894733466834395585465, -0.76585683089040890170196007602499,  0.44721359549995787219245357846376],
        [-0.52171527332952805089405501348665, -0.35779069914848704581089577914099,  0.07293249574894722364604149333900,  0.62807630128939140323751644245931,  0.44721359549995787219245357846376],
        [ 0.80478766909634613657686941223801, -0.35779069914848704581089577914099,  0.07293249574894722364604149333900,  0.13778052960101749846444363356568,  0.44721359549995787219245357846376],
      ],

      max_ulps = 1
    );

    Ok(())
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
    pretty_assert_ulps_eq!(
      r,
      array![
        [
          0.47098793460940813,
          0.778179233288187,
          0.7280398704422921,
          0.24949100191426288
        ],
        [
          0.4117271206212072,
          0.7624596970211635,
          0.7920856030429395,
          0.68138502792473
        ],
        [
          0.6713713388957947,
          0.8515178964314147,
          0.3118195942352078,
          0.7398254113552798
        ]
      ]
    );

    let r: Array2<f64> = random((3, 4), &mut rng);
    pretty_assert_ulps_eq!(
      r,
      array![
        [
          0.09054507828914282,
          0.7288312058240056,
          0.04697932861607845,
          0.5846614636962431
        ],
        [
          0.605341701141912,
          0.2310179777619814,
          0.667102895045365,
          0.9282527701718191
        ],
        [
          0.0715580637248685,
          0.2101276039082931,
          0.6010577625534657,
          0.94540503820318
        ]
      ]
    );
  }

  #[test]
  fn test_product_axis_axis_0() {
    let input = array![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]];
    let expected = array![15.0, 48.0];
    let result = product_axis(&input, Axis(0));
    pretty_assert_ulps_eq!(result, expected, epsilon = 1e-8);
  }

  #[test]
  fn test_product_axis_axis_1() {
    let input = array![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]];
    let expected = array![2.0, 12.0, 30.0];
    let result = product_axis(&input, Axis(1));
    pretty_assert_ulps_eq!(result, expected, epsilon = 1e-8);
  }
}
