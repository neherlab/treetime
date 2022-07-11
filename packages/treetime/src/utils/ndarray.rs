use eyre::Report;
use ndarray::{Array, Array1, Array2, Axis, Dimension, Ix2, NdProducer, RawData, RemoveAxis, ShapeBuilder, Zip};
use ndarray_rand::RandomExt;
use num_traits::real::Real;
use num_traits::{Bounded, Float, Num, NumCast};
use rand::distributions::uniform::SampleUniform;
use rand::distributions::Uniform;
use rand::Rng;
use std::ops::{Add, AddAssign, Deref, MulAssign};

pub fn to_col<T: Real>(a: &Array1<T>) -> Result<Array2<T>, Report> {
  Ok(a.to_shape((a.len(), 1))?.into_dimensionality::<Ix2>()?.to_owned())
}

pub fn to_row<T: Real>(b: &Array1<T>) -> Result<Array2<T>, Report> {
  Ok(b.to_shape((1, b.len()))?.into_dimensionality::<Ix2>()?.to_owned())
}

// Calculates outer product of 2 vectors
pub fn outer<T: 'static + Real>(a: &Array1<T>, b: &Array1<T>) -> Result<Array2<T>, Report> {
  let a = a.to_shape((a.len(), 1))?.into_dimensionality::<Ix2>()?;
  let b = b.to_shape((1, b.len()))?.into_dimensionality::<Ix2>()?;
  Ok(a.dot(&b))
}

/// Calculates min over given axis
#[inline]
pub fn min_axis<T: Real>(arr: &Array2<T>, axis: Axis) -> Result<Array1<T>, Report> {
  Ok(arr.fold_axis(axis, T::max_value(), |&a, &b| a.min(b)))
}

/// Calculates max over given axis
#[inline]
pub fn max_axis<T: Real>(arr: &Array2<T>, axis: Axis) -> Result<Array1<T>, Report> {
  Ok(arr.fold_axis(axis, T::min_value(), |&a, &b| a.max(b)))
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

/// Calculates cumulative sum over given axis
#[inline]
pub fn cumsum_axis<T: Copy + AddAssign, D: Dimension>(a: &Array<T, D>, axis: Axis) -> Array<T, D> {
  let mut result = a.to_owned();
  result.accumulate_axis_inplace(axis, |&prev, curr| *curr += prev);
  result
}

pub fn random<T: Copy + SampleUniform + NumCast, D: Dimension, Sh: ShapeBuilder<Dim = D>, R: Rng>(
  shape: Sh,
  rng: &mut R,
) -> Array<T, D> {
  let from: T = NumCast::from(0_i32).unwrap();
  let to: T = NumCast::from(1_i32).unwrap();
  Array::<T, D>::random_using(shape, Uniform::<T>::new::<T, T>(from, to), rng)
}

#[allow(clippy::excessive_precision, clippy::lossy_float_literal)]
#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use lazy_static::lazy_static;
  use ndarray::array;
  use ndarray_linalg::{Eigh, UPLO};
  use rand::SeedableRng;
  use rand_isaac::Isaac64Rng;
  use rstest::rstest;
  use crate::pretty_assert_eq;

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

    assert_ulps_eq!(cumsum_axis(&INPUT, Axis(0)), expected);
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

    assert_ulps_eq!(cumsum_axis(&INPUT, Axis(1)), expected);
  }

  #[rstest]
  fn computes_outer_product() -> Result<(), Report> {
    assert_ulps_eq!(
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
    //
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
    assert_ulps_eq!(
      eigvals,

      // Rust ndarray_linalg::eigh() result:
      //    [ -1.2500000000000002, -1.25,            -1.2499999999999998,  -1.249999999999999,   5.551115123125783e-17 ]

      // NumPy np.linalg.eigh() result:
      array![ -1.25000000e+00,     -1.25000000e+00,  -1.25000000e+00,      -1.25000000e+00,      5.55111512e-17        ]
    );

    #[rustfmt::skip]
    assert_ulps_eq!(
      eigvecs,

      // Rust ndarray_linalg::eigh(UPLO::Lower) result:
      // [
      //   [  0.0,                     -0.6779192194392384,    0.5834599659915785,  -0.0,                    -0.447213595499958   ],
      //   [ -4.3301795267950775e-17,  -0.39545287800622264,  -0.8022574532384199,  -1.1715368998076453e-16, -0.4472135954999581  ],
      //   [ -0.28307239576681803,      0.35779069914848716,   0.0729324957489474,  -0.7658568308904088,     -0.4472135954999579  ],
      //   [ -0.521715273329528,        0.357790699148487,     0.07293249574894722,  0.6280763012893915,     -0.44721359549995787 ],
      //   [  0.8047876690963461,       0.357790699148487,     0.07293249574894722,  0.1377805296010175,     -0.44721359549995787 ]
      // ]

      // NumPy np.linalg.eigh() result:
      array![
           [  0.00000000e+00,          -6.77919219e-01,        5.83459966e-01,      -0.00000000e+00,         -4.47213595e-01      ],
           [ -4.16333634e-17,          -3.95452878e-01,       -8.02257453e-01,       4.16333634e-17,         -4.47213595e-01      ],
           [ -2.83072396e-01,           3.57790699e-01,        7.29324957e-02,      -7.65856831e-01,         -4.47213595e-01      ],
           [ -5.21715273e-01,           3.57790699e-01,        7.29324957e-02,       6.28076301e-01,         -4.47213595e-01      ],
           [  8.04787669e-01,           3.57790699e-01,        7.29324957e-02,       1.37780530e-01,         -4.47213595e-01      ],
      ]
    );

    Ok(())
  }

  #[rstest]
  fn generates_predictable_random_uniform() {
    let mut rng = Isaac64Rng::seed_from_u64(42);

    let r: Array2<f64> = random((3, 4), &mut rng);
    assert_ulps_eq!(
      r,
      array![
        [0.95799470, 0.47098792, 0.39911770, 0.77817917],
        [0.73372230, 0.72803986, 0.86803260, 0.24949098],
        [0.94433236, 0.41172707, 0.22396588, 0.76245964],
      ]
    );

    let r: Array2<f64> = random((3, 4), &mut rng);
    assert_ulps_eq!(
      r,
      array![
        [0.70186340, 0.79208550, 0.92076266, 0.68138490],
        [0.02917111, 0.67137120, 0.13164222, 0.85151780],
        [0.49220824, 0.31181955, 0.52067375, 0.73982537],
      ]
    );
  }
}
