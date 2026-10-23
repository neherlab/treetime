use crate::gtr::gtr::GTR;
use eyre::Report;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum InferredRate {
  Fitted,
  Unit,
}

impl InferredRate {
  pub(crate) const fn for_partition_count(partition_count: usize) -> Self {
    if partition_count == 1 { Self::Unit } else { Self::Fitted }
  }

  pub(crate) fn apply(self, gtr: &mut GTR) -> Result<(), Report> {
    match self {
      Self::Fitted => Ok(()),
      Self::Unit => gtr.set_mu(1.0),
    }
  }
}
