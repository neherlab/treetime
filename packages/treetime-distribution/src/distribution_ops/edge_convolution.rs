use crate::Distribution;
use crate::DistributionFunction;
use crate::distribution_ops::convolve::distribution_convolution_fine;
use crate::distribution_ops::mass_domain::{
  SpacingLimit, mass_bounded_domain, mass_sizable_profile, peak_normalized_if_mass_sizable, resample_to_mass_window,
};
use crate::distribution_ops::mass_profile::TailMass;
use crate::policy::NegLog;
use eyre::Report;
use ordered_float::OrderedFloat;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, MaxGridPoints, Side};
use treetime_utils::cast::ToF64Exact;

pub fn convolve_across_edge(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
  soft: Side,
  eps: f64,
  grid_points: usize,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  let conv = WindowedConvolution::new(a, b, soft, eps, max_points)?;
  match soft {
    Side::Right => conv.resample(grid_points, max_points, OutputSpacing::Capped),
    Side::Left => conv.resample(grid_points, max_points, OutputSpacing::Fine),
  }
}

pub fn convolve_forward_across_edge(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
  eps: f64,
  grid_points: usize,
  max_points: MaxGridPoints,
) -> Result<ForwardConvolution, Report> {
  let conv = WindowedConvolution::new(a, b, Side::Right, eps, max_points)?;
  Ok(ForwardConvolution {
    conv,
    grid_points,
    max_points,
  })
}

pub struct ForwardConvolution {
  conv: WindowedConvolution,
  grid_points: usize,
  max_points: MaxGridPoints,
}

impl ForwardConvolution {
  pub fn capped(&self) -> Result<Distribution<NegLog>, Report> {
    self
      .conv
      .resample(self.grid_points, self.max_points, OutputSpacing::Capped)
  }

  pub fn resampled_to_resolve(
    &self,
    posterior: &Distribution<NegLog>,
    eps: f64,
  ) -> Result<Option<Distribution<NegLog>>, Report> {
    let (
      WindowedConvolution::Window {
        normalized,
        window: (lo, hi),
        cap: Some(cap),
      },
      Distribution::Function(f),
    ) = (&self.conv, posterior)
    else {
      return Ok(None);
    };
    let Some(profile) = mass_sizable_profile(f) else {
      return Ok(None);
    };
    let Some(fractions) = profile.fractions(TailMass::Required) else {
      return Ok(None);
    };
    let intervals = self.grid_points.saturating_sub(1).max(1).to_f64_exact();
    let capped_dx = ((hi - lo) / intervals).min(*cap);
    let mass_width = fractions.survival_inverse(eps) - fractions.quantile(eps);
    let resolving_dx = normalized.dx().max(mass_width / intervals);
    if capped_dx <= resolving_dx {
      return Ok(None);
    }
    resample_to_mass_window(
      normalized,
      *lo,
      *hi,
      self.grid_points,
      self.max_points,
      SpacingLimit::AtMost(resolving_dx),
    )
    .map(Some)
  }
}

enum WindowedConvolution {
  Exact(Distribution<NegLog>),
  Window {
    normalized: DistributionFunction<f64, NegLog>,
    window: (f64, f64),
    cap: Option<f64>,
  },
}

impl WindowedConvolution {
  fn new(
    a: &Distribution<NegLog>,
    b: &Distribution<NegLog>,
    soft: Side,
    eps: f64,
    max_points: MaxGridPoints,
  ) -> Result<Self, Report> {
    let is_fft = matches!((a, b), (Distribution::Function(_), Distribution::Function(_)));
    let conv = distribution_convolution_fine(a, b, max_points)?;
    let conv = if is_fft || has_no_soft_law(&conv, soft) {
      conv.fit_soft_tail(soft, DEFAULT_TAIL_FIT_POINTS)?
    } else {
      conv
    };
    let conv = match soft {
      Side::Left => conv.with_right_extrap(BoundaryBehavior::Hard)?,
      Side::Right => conv.with_left_extrap(BoundaryBehavior::Hard)?,
    };
    let Distribution::Function(conv) = conv else {
      return Ok(Self::Exact(conv));
    };
    if !is_fft {
      return Ok(Self::Exact(Distribution::Function(conv)));
    }
    let Some(sized) = peak_normalized_if_mass_sizable(&conv) else {
      return Ok(Self::Exact(Distribution::Function(conv)));
    };

    let window = match convolution_output_window(a, b, eps)? {
      Some((mut lo, mut hi)) => {
        match soft {
          Side::Left => hi = hi.min(sized.normalized.x_max()),
          Side::Right => lo = lo.max(sized.normalized.x_min()),
        }
        if hi > lo { (lo, hi) } else { sized.bounded_domain(eps)? }
      },
      None => sized.bounded_domain(eps)?,
    };

    Ok(Self::Window {
      normalized: sized.normalized,
      window,
      cap: coarser_operand_spacing(a, b),
    })
  }

  fn resample(
    &self,
    grid_points: usize,
    max_points: MaxGridPoints,
    spacing: OutputSpacing,
  ) -> Result<Distribution<NegLog>, Report> {
    match self {
      Self::Exact(conv) => conv.repeak(),
      Self::Window {
        normalized,
        window: (lo, hi),
        cap,
      } => {
        let spacing = match (spacing, cap) {
          (OutputSpacing::Capped, Some(max_dx)) => SpacingLimit::AtMost(*max_dx),
          (OutputSpacing::Capped, None) | (OutputSpacing::Fine, _) => SpacingLimit::Input,
        };
        resample_to_mass_window(normalized, *lo, *hi, grid_points, max_points, spacing)
      },
    }
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum OutputSpacing {
  Capped,
  Fine,
}

fn coarser_operand_spacing(a: &Distribution<NegLog>, b: &Distribution<NegLog>) -> Option<f64> {
  [a, b]
    .into_iter()
    .filter_map(|dist| match dist {
      Distribution::Function(f) => Some(f.dx()),
      Distribution::Empty | Distribution::Point(_) | Distribution::Range(_) | Distribution::Formula(_) => None,
    })
    .max_by_key(|dx| OrderedFloat(*dx))
}

fn has_no_soft_law(conv: &Distribution<NegLog>, soft: Side) -> bool {
  let Distribution::Function(f) = conv else {
    return false;
  };
  let law = match soft {
    Side::Left => f.left_extrap(),
    Side::Right => f.right_extrap(),
  };
  law == BoundaryBehavior::Error
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
