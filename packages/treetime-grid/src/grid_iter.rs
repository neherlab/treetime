use crate::InterpElem;
use crate::grid::Grid;
use num_traits::Float;

/// Iterator over grid x coordinates
pub struct GridIter<T: InterpElem> {
  grid: Grid<T>,
  index: usize,
}

impl<T: InterpElem> GridIter<T> {
  pub(crate) fn new(grid: Grid<T>) -> Self {
    Self { grid, index: 0 }
  }
}

impl<T: InterpElem> Iterator for GridIter<T>
where
  T: Float,
{
  type Item = T;

  fn next(&mut self) -> Option<Self::Item> {
    (self.index < self.grid.len()).then(|| {
      let x = self.grid.x_at(self.index);
      self.index += 1;
      x
    })
  }

  fn size_hint(&self) -> (usize, Option<usize>) {
    let remaining = self.grid.len() - self.index;
    (remaining, Some(remaining))
  }
}

impl<T: InterpElem> ExactSizeIterator for GridIter<T>
where
  T: Float,
{
  fn len(&self) -> usize {
    self.grid.len() - self.index
  }
}
