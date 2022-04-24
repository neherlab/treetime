use crate::make_report;
use crate::utils::ndarray::max_axis;
use eyre::Report;
use ndarray::prelude::*;
use ndarray_einsum_beta::einsum;
use num_traits::real::Real;

/// Returns a normalized version of a profile matrix
pub fn normalize_profile(in_prof: &Array2<f32>, log: bool) -> Result<Array2<f32>, Report> {
  let tmp_prof = if log {
    let tmp_prefactor = max_axis(in_prof, Axis(1))?;

    let mut tmp_prof = in_prof.t().to_owned() - tmp_prefactor;
    tmp_prof = tmp_prof.mapv(f32::exp);

    tmp_prof
  } else {
    in_prof.clone()
  };

  let norm_vector = tmp_prof.sum_axis(Axis(1));
  let one_div_norm_vector = 1.0 / norm_vector;

  let norm_prof = einsum("ai,a->ai", &[&tmp_prof, &one_div_norm_vector])
    .map_err(|err| make_report!("einsum: {err}"))?
    .into_dimensionality::<Ix2>()?;

  Ok(norm_prof)
}
