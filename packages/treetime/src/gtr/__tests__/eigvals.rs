#[cfg(test)]
mod tests {
  use crate::gtr::gtr::GTR;
  use ndarray::Array1;

  impl GTR {
    pub(crate) fn eigvals(&self) -> &Array1<f64> {
      &self.eigvals
    }
  }
}
