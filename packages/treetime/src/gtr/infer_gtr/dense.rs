use crate::gtr::gtr::GTR;
use crate::make_error;
use eyre::Report;

/// Infer GTR model from dense partition data.
///
/// Dense representation uses full probability matrices. GTR inference from dense
/// requires computing expected mutation counts from probability distributions.
/// Not yet implemented.
pub fn infer_gtr_dense() -> Result<GTR, Report> {
  make_error!("GTR model inference is not yet implemented for dense representation")
}
