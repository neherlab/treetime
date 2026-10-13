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

#[expect(
  clippy::unwrap_used,
  reason = "documented non-empty precondition; into_option is Some for a non-empty iterator"
)]
pub fn minmax<'a, T: PartialOrd + Clone + 'a, I: Iterator<Item = &'a T> + 'a>(iter: I) -> (T, T) {
  let (x_min, x_max) = iter.minmax().into_option().unwrap();
  (x_min.clone(), x_max.clone())
}
