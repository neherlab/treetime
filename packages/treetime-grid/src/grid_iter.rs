use crate::grid::Grid;

pub struct GridIter {
  grid: Grid,
  index: usize,
}

impl GridIter {
  pub(crate) fn new(grid: Grid) -> Self {
    Self { grid, index: 0 }
  }
}

impl Iterator for GridIter {
  type Item = f64;

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
