use approx::{UlpsEq, ulps_eq};
use eyre::Report;
use itertools::Itertools;
use ndarray::{
  Array, Array1, Array2, ArrayBase, ArrayView, Axis, Data, Dimension, Ix1, Ix2, RemoveAxis, ShapeBuilder, ShapeError,
  Zip, s, stack,
};
use ndarray_rand::RandomExt;
use ndarray_rand::rand::Rng;
use ndarray_rand::rand::distributions::Uniform;
use ndarray_rand::rand::distributions::uniform::SampleUniform;
use num_traits::real::Real;
use num_traits::{Bounded, Float, NumCast, One, Zero};
use std::f64::consts::E;
use std::ops::{AddAssign, Mul};

pub fn first<T: Copy>(a: &Array1<T>) -> T {
  a[0]
}

pub fn last<T: Copy>(a: &Array1<T>) -> T {
  a[a.len() - 1]
}

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

/// Pad array to target length with zeros on the right, or truncate if input is longer
pub fn ndarray_pad_zeros_right(input: &Array1<f64>, target_length: usize) -> Array1<f64> {
  let mut result = Array1::zeros(target_length);
  let copy_len = input.len().min(target_length);
  result.slice_mut(s![..copy_len]).assign(&input.slice(s![..copy_len]));
  result
}

/// Create uniform grid with specified start, spacing, and length
pub fn ndarray_uniform_grid(start: f64, spacing: f64, length: usize) -> Array1<f64> {
  Array1::linspace(start, start + spacing * ((length - 1) as f64), length)
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

/// Calculate product over given axis.
/// Backported from ndarray 0.16.1.
/// https://github.com/rust-ndarray/ndarray/blob/6f77377d7d508550bf516e54c142cee3ab243aeb/src/numeric/impl_numeric.rs#L267-L282
pub fn product_axis<A, D>(a: &Array<A, D>, axis: Axis) -> Array<A, D::Smaller>
where
  A: Clone + One + Mul<Output = A>,
  D: RemoveAxis,
{
  let min_stride_axis = a.raw_dim().min_stride_axis(&a.raw_dim());
  if axis == min_stride_axis {
    Zip::from(a.lanes(axis)).map_collect(|lane| lane.product())
  } else {
    let mut res = Array::ones(a.raw_dim().remove_axis(axis));
    for subview in a.axis_iter(axis) {
      res = res * &subview;
    }
    res
  }
}

pub fn random<T: Copy + SampleUniform + NumCast, D: Dimension, Sh: ShapeBuilder<Dim = D>, R: Rng>(
  shape: Sh,
  rng: &mut R,
) -> Array<T, D> {
  let from: T = NumCast::from(0_i32).unwrap();
  let to: T = NumCast::from(1_i32).unwrap();
  Array::<T, D>::random_using(shape, Uniform::<T>::new::<T, T>(from, to), rng)
}

/// Reverse 1D array in place by inverting axis 0
pub fn reverse_inplace<T>(arr: &mut Array1<T>) {
  arr.invert_axis(Axis(0));
}

/// Reverse 1D array by inverting axis 0
pub fn reverse<T: Clone, S: Data<Elem = T>>(arr: &ArrayBase<S, Ix1>) -> Array1<T> {
  let mut reversed = arr.view().to_owned();
  reverse_inplace(&mut reversed);
  reversed
}

/// Check if array spacing is uniform
pub fn has_uniform_spacing<T: Float + UlpsEq>(grid: &Array1<T>) -> bool {
  if grid.len() < 2 {
    return true;
  }

  let spacing = grid[1] - grid[0];
  grid
    .windows(2)
    .into_iter()
    .all(|w| ulps_eq!(w[1] - w[0], spacing, max_ulps = 100))
}
