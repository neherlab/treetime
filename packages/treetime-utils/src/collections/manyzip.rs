/// Zip an unknown number of iterators
///
/// Borrowed with modifications from https://stackoverflow.com/a/55292215
/// Author: Shepmaster
pub struct Manyzip<T>(pub Vec<T>);

impl<T> Iterator for Manyzip<T>
where
  T: Iterator,
{
  type Item = Vec<T::Item>;

  fn next(&mut self) -> Option<Self::Item> {
    self.0.iter_mut().map(Iterator::next).collect()
  }
}
