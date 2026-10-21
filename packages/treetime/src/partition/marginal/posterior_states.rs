use crate::constants::POSTERIOR_TIE_TOLERANCE;
use ndarray::{ArrayBase, Data, Ix1};
use treetime_utils::array::ndarray::argmax_first_within;

pub fn most_likely_state<S: Data<Elem = f64>>(posterior: &ArrayBase<S, Ix1>) -> Option<usize> {
  argmax_first_within(posterior, POSTERIOR_TIE_TOLERANCE)
}
