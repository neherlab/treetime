use crate::make_error;
use eyre::Report;
use itertools::Itertools;

pub fn get_one<T>(x: &[T]) -> Option<&T> {
  x.first()
}

pub fn get_exactly_one<T>(x: &[T]) -> Result<&T, Report> {
  match x.len() {
    1 => Ok(&x[0]),
    _ => make_error!("Expected exactly one element, but found {}", x.len()),
  }
}

pub fn get_exactly_one_mut<T>(x: &mut [T]) -> Result<&mut T, Report> {
  match x.len() {
    1 => Ok(&mut x[0]),
    _ => make_error!("Expected exactly one element, but found {}", x.len()),
  }
}

/// Find min and max of a non-empty iterator
pub fn minmax<'a, T: PartialOrd + Clone + 'a, I: Iterator<Item = &'a T> + 'a>(iter: I) -> (T, T) {
  let (x_min, x_max) = iter.minmax().into_option().unwrap();
  (x_min.clone(), x_max.clone())
}
