use eyre::Report;
use ndarray::{Array, Array1, Array2, Axis, Dim, Dimension, Ix1, Ix2, NdProducer, RawData, Zip};
use num_traits::Bounded;
use std::ops::MulAssign;

pub fn to_col(a: &Array1<f32>) -> Result<Array2<f32>, Report> {
  Ok(a.to_shape((a.len(), 1))?.into_dimensionality::<Ix2>()?.to_owned())
}

pub fn to_row(b: &Array1<f32>) -> Result<Array2<f32>, Report> {
  Ok(b.to_shape((1, b.len()))?.into_dimensionality::<Ix2>()?.to_owned())
}

// Calculates outer product of 2 vectors
pub fn outer(a: &Array1<f32>, b: &Array1<f32>) -> Result<Array2<f32>, Report> {
  let a = a.to_shape((a.len(), 1))?.into_dimensionality::<Ix2>()?;
  let b = b.to_shape((1, b.len()))?.into_dimensionality::<Ix2>()?;
  Ok(a.dot(&b))
}

/// Calculates min sum over given axis
#[inline]
pub fn min_axis(arr: &Array2<f32>, axis: Axis) -> Result<Array1<f32>, Report> {
  Ok(arr.fold_axis(axis, f32::MAX, |a, b| a.min(*b)))
}

/// Calculates max sum over given axis
#[inline]
pub fn max_axis(arr: &Array2<f32>, axis: Axis) -> Result<Array1<f32>, Report> {
  Ok(arr.fold_axis(axis, f32::MIN, |a, b| a.max(*b)))
}

/// Finds index of min value over given axis
#[inline]
pub fn argmin_axis<T: Copy + Ord + PartialOrd + Bounded>(arr: &Array1<T>, axis: Axis) -> usize {
  let mut i_curr = 0;
  let res = arr
    .fold_axis(axis, (0_usize, T::max_value()), |(i_min, x_min), x| {
      let res = if x < x_min { (i_curr, *x) } else { (*i_min, *x_min) };
      i_curr += 1;
      res
    })
    .into_raw_vec();
  res[0].0
}

/// Finds index of max value over given axis
#[inline]
pub fn argmax_axis<T: Copy + Ord + PartialOrd + Bounded>(arr: &Array1<T>, axis: Axis) -> usize {
  let mut i_curr = 0;
  let res = arr
    .fold_axis(axis, (0_usize, T::min_value()), |(i_max, x_max), x| {
      let res = if x > x_max { (i_curr, *x) } else { (*i_max, *x_max) };
      i_curr += 1;
      res
    })
    .into_raw_vec();
  res[0].0
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
pub fn cumsum_axis(profile: &Array2<f32>, axis: Axis) -> Result<Array1<f32>, Report> {
  let mut result = profile.to_owned();
  result.accumulate_axis_inplace(axis, |&prev, curr| *curr += prev);
  let result = result.into_dimensionality::<Ix1>()?;
  Ok(result)
}

#[allow(clippy::excessive_precision)]
#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;
  use ndarray_linalg::{Eigh, UPLO};
  use rstest::rstest;

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
}
