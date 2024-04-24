use crate::make_internal_error;
use eyre::Report;

#[macro_export]
macro_rules! vec_of_owned {
    () => (
        vec![]
    );
    ($($x:expr),+ $(,)?) => (
        vec![$($x),+].into_iter().map(|x| x.to_owned()).collect_vec()
    );
}

pub fn get_exactly_one<T>(arr: &[T]) -> Result<&T, Report> {
  if arr.len() != 1 {
    make_internal_error!("Expected exactly one element, but found '{}'", arr.len())
  } else {
    Ok(&arr[0])
  }
}
