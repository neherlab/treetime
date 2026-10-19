use crate::Distribution;
use crate::distribution_ops::convolve::distribution_convolution_fine;
use crate::distribution_ops::mass_domain::{
  experiment_coarse, experiment_trim_mode, mass_bounded_domain, mass_sizable_profile, mass_trimmed_domain,
  peak_normalized_if_mass_sizable, resample_to_mass_window, resample_trimmed,
};
use crate::policy::NegLog;
use eyre::Report;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, MaxGridPoints, Side};

pub fn convolve_across_edge(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
  soft: Side,
  eps: f64,
  grid_points: usize,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  let (a, b) = if std::env::var_os("TT_TAILS").is_some() {
    (&extend_soft_tails(a)?, &extend_soft_tails(b)?)
  } else {
    (a, b)
  };
  let conv = distribution_convolution_fine(a, b, max_points)?.fit_soft_tail(soft, DEFAULT_TAIL_FIT_POINTS)?;
  let conv = match soft {
    Side::Left => conv.with_right_extrap(BoundaryBehavior::Hard)?,
    Side::Right => conv.with_left_extrap(BoundaryBehavior::Hard)?,
  };
  let Distribution::Function(conv) = conv else {
    return conv.repeak();
  };
  let Some(sized) = peak_normalized_if_mass_sizable(&conv) else {
    return Distribution::Function(conv).repeak();
  };
  let normalized = &sized.normalized;

  if let Some(mode) = experiment_trim_mode().filter(|_| soft == Side::Right) {
    // `p`: the posterior alone is trimmed; the convolution keeps its hard side and only takes the coarser spacing
    let operand_domain = if mode == b'p' { operand_mass_domain } else { trimmed_operand_domain };
    let window = match (operand_domain(a, eps)?, operand_domain(b, eps)?) {
      (Some((lo_a, hi_a)), Some((lo_b, hi_b))) => Some((lo_a + lo_b, hi_a + hi_b)),
      _ => None,
    };
    let (lo, hi) = match window {
      Some((lo, hi)) => {
        let (lo, hi) = (lo.max(normalized.x_min()), hi.min(normalized.x_max()));
        if hi > lo {
          (lo, hi)
        } else if mode == b'p' {
          sized.bounded_domain(eps)?
        } else {
          mass_trimmed_domain(&sized.profile, eps)?
        }
      },
      None if mode == b'p' => sized.bounded_domain(eps)?,
      None => mass_trimmed_domain(&sized.profile, eps)?,
    };
    let max_dx = if experiment_coarse() {
      match (a, b) {
        (Distribution::Function(fa), Distribution::Function(fb)) => Some(fa.dx().max(fb.dx())),
        (Distribution::Function(f), _) | (_, Distribution::Function(f)) => Some(f.dx()),
        _ => None,
      }
    } else {
      None
    };
    return resample_trimmed(normalized, lo, hi, grid_points, max_points, max_dx, mode);
  }

  let (lo, hi) = match convolution_output_window(a, b, eps)? {
    Some((mut lo, mut hi)) => {
      match soft {
        Side::Left => hi = hi.min(normalized.x_max()),
        Side::Right => lo = lo.max(normalized.x_min()),
      }
      if hi > lo { (lo, hi) } else { sized.bounded_domain(eps)? }
    },
    None => sized.bounded_domain(eps)?,
  };

  resample_to_mass_window(normalized, lo, hi, grid_points, max_points)
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
    Distribution::Function(f) => mass_sizable_profile(f)
      .map(|profile| mass_bounded_domain(&profile, eps))
      .transpose(),
    Distribution::Empty | Distribution::Formula(_) => Ok(None),
  }
}

fn trimmed_operand_domain(dist: &Distribution<NegLog>, eps: f64) -> Result<Option<(f64, f64)>, Report> {
  match dist {
    Distribution::Function(f) => mass_sizable_profile(f)
      .map(|profile| mass_trimmed_domain(&profile, eps))
      .transpose(),
    other => operand_mass_domain(other, eps),
  }
}

/// EXPERIMENT (kb issue M-distribution-convolution-drops-operand-tails): puts the soft (`Linear`) tails of a function
/// onto its grid for `TAIL_DECAY_NATS` of decay, at the grid's own spacing, so that the convolution sees their mass
fn extend_soft_tails(dist: &Distribution<NegLog>) -> Result<Distribution<NegLog>, Report> {
  const TAIL_DECAY_NATS: f64 = 12.0;
  let Distribution::Function(f) = dist else {
    return Ok(dist.clone());
  };
  let reach = |b: BoundaryBehavior| match b {
    BoundaryBehavior::Linear(law) if law.slope.abs() > 0.0 => TAIL_DECAY_NATS / law.slope.abs(),
    _ => 0.0,
  };
  let (lo, hi) = (f.x_min() - reach(f.left_extrap()), f.x_max() + reach(f.right_extrap()));
  if lo == f.x_min() && hi == f.x_max() {
    return Ok(dist.clone());
  }
  let n_points = ((hi - lo) / f.dx()).ceil() as usize + 1;
  Ok(Distribution::Function(f.resample_range_n_points((lo, hi), n_points)?))
}
