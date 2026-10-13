use crate::cast::ToF64Exact;
use eyre::Report;
use ndarray::{Array, Array1, Array2, ArrayBase, ArrayView, Axis, Data, Dimension, Ix1, Ix2, ShapeBuilder, Zip};
use ndarray_rand::RandomExt;
use ndarray_rand::rand::Rng;
use ndarray_rand::rand::distributions::Uniform;
use ndarray_rand::rand::distributions::uniform::SampleUniform;
use ndarray_stats::QuantileExt;
use num_traits::float::FloatCore;
use num_traits::real::Real;
use num_traits::{Float, One, Zero};
use ordered_float::OrderedFloat;
use std::f64::consts::E;
use std::ops::AddAssign;

pub fn first<T: Copy>(a: &Array1<T>) -> T {
  a[0]
}

pub fn last<T: Copy>(a: &Array1<T>) -> T {
  a[a.len() - 1]
}

pub fn outer<T: 'static + Real>(a: &Array1<T>, b: &Array1<T>) -> Result<Array2<T>, Report> {
  let a = a.to_shape((a.len(), 1))?.into_dimensionality::<Ix2>()?;
  let b = b.to_shape((1, b.len()))?.into_dimensionality::<Ix2>()?;
  Ok(a.dot(&b))
}

pub fn argmax_first<T: PartialOrd>(arr: &ArrayView<T, Ix1>) -> Option<usize> {
  if arr.is_empty() {
    return None;
  }
  let mut max_idx = 0;
  let mut max_val = &arr[0];
  for (i, val) in arr.iter().enumerate().skip(1) {
    if val > max_val {
      max_idx = i;
      max_val = val;
    }
  }
  Some(max_idx)
}

#[inline]
pub fn max_or<S: Data<Elem = f64>>(arr: &ArrayBase<S, Ix1>, default: f64) -> f64 {
  arr.max().map_or(default, |&max| max)
}

#[inline]
pub fn min_or<S: Data<Elem = f64>>(arr: &ArrayBase<S, Ix1>, default: f64) -> f64 {
  arr.min().map_or(default, |&min| min)
}

#[inline]
pub fn is_max_above<S: Data<Elem = f64>>(arr: &ArrayBase<S, Ix1>, threshold: f64) -> bool {
  arr.max().is_ok_and(|&max| max >= threshold)
}

pub fn minimum<T: Copy + PartialOrd, D: Dimension>(x: &Array<T, D>, y: &Array<T, D>) -> Array<T, D> {
  assert_eq!(x.shape(), y.shape());
  Zip::from(x).and(y).map_collect(|&a, &b| if a < b { a } else { b })
}

pub fn maximum<T: Copy + PartialOrd, D: Dimension>(x: &Array<T, D>, y: &Array<T, D>) -> Array<T, D> {
  assert_eq!(x.shape(), y.shape());
  Zip::from(x).and(y).map_collect(|&a, &b| if a > b { a } else { b })
}

pub fn clamp_min<T: Copy + PartialOrd, D: Dimension>(a: &Array<T, D>, lower: T) -> Array<T, D> {
  a.mapv(|x| num_traits::clamp_min(x, lower))
}

pub fn clamp_max<T: Copy + PartialOrd, D: Dimension>(a: &Array<T, D>, upper: T) -> Array<T, D> {
  a.mapv(|x| num_traits::clamp_max(x, upper))
}

pub fn clamp<T: Copy + PartialOrd, D: Dimension>(a: &Array<T, D>, lower: T, upper: T) -> Array<T, D> {
  a.mapv(|x| num_traits::clamp(x, lower, upper))
}

pub fn exp<T, S, D>(arr: &ArrayBase<S, D>) -> Array<T, D>
where
  T: Float + From<f64>,
  S: Data<Elem = T>,
  D: Dimension,
{
  arr.mapv(|x| Float::exp(x))
}

pub fn log<T: Float + From<f64>, D: Dimension>(arr: &Array<T, D>) -> Array<T, D> {
  arr.mapv(|x| Float::log(x, E.into()))
}

pub fn zeros<T, D, Sh>(shape: Sh) -> Array<T, D>
where
  T: Zero + Clone,
  D: Dimension,
  Sh: ShapeBuilder<Dim = D>,
{
  Array::<T, D>::zeros(shape)
}

#[inline]
pub fn cumsum_axis<T: Copy + AddAssign, D: Dimension>(a: &Array<T, D>, axis: Axis) -> Array<T, D> {
  let mut result = a.to_owned();
  result.accumulate_axis_inplace(axis, |&prev, curr| *curr += prev);
  result
}

pub fn random<T: Copy + SampleUniform + Zero + One, D: Dimension, Sh: ShapeBuilder<Dim = D>, R: Rng>(
  shape: Sh,
  rng: &mut R,
) -> Array<T, D> {
  Array::<T, D>::random_using(shape, Uniform::<T>::new::<T, T>(T::zero(), T::one()), rng)
}

pub fn reverse<T: Clone, S: Data<Elem = T>>(arr: &ArrayBase<S, Ix1>) -> Array1<T> {
  let mut reversed = arr.view().to_owned();
  reverse_inplace(&mut reversed);
  reversed
}

pub fn reverse_inplace<T>(arr: &mut Array1<T>) {
  arr.invert_axis(Axis(0));
}

pub fn sorted<T: FloatCore + Clone>(arr: &Array1<T>) -> Array1<T> {
  let mut result = arr.clone();
  sort_inplace(&mut result);
  result
}

pub fn partition_point<T>(arr: &Array1<T>, mut pred: impl FnMut(&T) -> bool) -> usize {
  if let Some(values) = arr.as_slice() {
    values.partition_point(pred)
  } else {
    arr.iter().take_while(|&value| pred(value)).count()
  }
}

pub fn sort_inplace<T: FloatCore>(arr: &mut Array1<T>) {
  if let Some(values) = arr.as_slice_mut() {
    values.sort_by_key(|x| OrderedFloat(*x));
  } else {
    *arr = Array1::from_iter(arr.iter().copied().sorted_by_key(|x| OrderedFloat(*x)));
  }
}

const MAX_SPACING_ULPS: f64 = 64.0;

pub fn has_uniform_spacing<T: Float + From<f64>>(grid: &Array1<T>) -> bool {
  if grid.len() < 2 {
    return true;
  }

  let spacing = grid[1] - grid[0];
  let scale = grid[0].abs().max(grid[grid.len() - 1].abs());
  let tolerance = <T as From<f64>>::from(MAX_SPACING_ULPS) * scale * T::epsilon();
  grid
    .windows(2)
    .into_iter()
    .all(|w| (w[1] - w[0] - spacing).abs() <= tolerance)
}
