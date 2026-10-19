use crate::gtr::gtr::GTR;
use eyre::Report;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum InferredRate {
  Fitted,
  Unit,
}

impl InferredRate {
  pub(crate) const fn for_file_count(file_count: usize) -> Self {
    if file_count == 1 { Self::Unit } else { Self::Fitted }
  }

  pub(crate) const fn for_separate_files() -> Self {
    Self::for_file_count(1)
  }

  pub(crate) fn apply(self, gtr: &mut GTR) -> Result<(), Report> {
    match self {
      Self::Fitted => Ok(()),
      Self::Unit => gtr.set_mu(1.0),
    }
  }
}
