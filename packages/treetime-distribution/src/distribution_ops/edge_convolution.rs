use crate::Distribution;
use crate::distribution_ops::convolve::distribution_convolution_fine;
use crate::distribution_ops::mass_domain::{
  mass_bounded_domain, peak_normalized_if_mass_sizable, resample_to_mass_window,
};
use crate::policy::NegLog;
use eyre::Report;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, Side};

pub fn convolve_across_edge(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
  soft: Side,
  eps: f64,
  grid_points: usize,
) -> Result<Distribution<NegLog>, Report> {
  let conv = distribution_convolution_fine(a, b)?.fit_soft_tail(soft, DEFAULT_TAIL_FIT_POINTS)?;
  let conv = match soft {
    Side::Left => conv.with_right_extrap(BoundaryBehavior::Hard)?,
    Side::Right => conv.with_left_extrap(BoundaryBehavior::Hard)?,
  };
  let Distribution::Function(conv) = conv else {
    return conv.normalize();
  };
  let Some(normalized) = peak_normalized_if_mass_sizable(&conv) else {
    return Distribution::Function(conv).normalize();
  };

  let (lo, hi) = match convolution_output_window(a, b, eps)? {
    Some((mut lo, mut hi)) => {
      match soft {
        Side::Left => hi = hi.min(normalized.x_max()),
        Side::Right => lo = lo.max(normalized.x_min()),
      }
      if hi > lo {
        (lo, hi)
      } else {
        mass_bounded_domain(&normalized, eps)?
      }
    },
    None => mass_bounded_domain(&normalized, eps)?,
  };

  resample_to_mass_window(&normalized, lo, hi, grid_points)
}

fn convolution_output_window(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
  eps: f64,
) -> Result<Option<(f64, f64)>, Report> {
  let (Some((lo_a, hi_a)), Some((lo_b, hi_b))) = (operand_mass_domain(a, eps)?, operand_mass_domain(b, eps)?) else {
    return Ok(None);
  };
  Ok(Some((lo_a + lo_b, hi_a + hi_b)))
}

fn operand_mass_domain(dist: &Distribution<NegLog>, eps: f64) -> Result<Option<(f64, f64)>, Report> {
  match dist {
    Distribution::Point(p) => Ok(Some((p.t(), p.t()))),
    Distribution::Range(r) => Ok(Some((r.start(), r.end()))),
    Distribution::Function(f) => peak_normalized_if_mass_sizable(f)
      .map(|normalized| mass_bounded_domain(&normalized, eps))
      .transpose(),
    Distribution::Empty | Distribution::Formula(_) => Ok(None),
  }
}
