use crate::gtr::gtr::GTR;
use ndarray::{Array2, Axis};

pub(super) fn rate_matrix(gtr: &GTR) -> Array2<f64> {
  let mut q = (&gtr.W * &gtr.pi).t().to_owned();
  let diag = -q.sum_axis(Axis(0));
  q.diag_mut().assign(&diag);
  q
}
