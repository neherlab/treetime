use crate::Distribution;
use crate::DistributionFunction;
use crate::distribution_ops::convolve::distribution_convolution_fine;
use crate::distribution_ops::mass_domain::{
  SpacingLimit, mass_bounded_domain, mass_sizable_profile, peak_normalized_if_mass_sizable, resample_to_mass_window,
  trim_to_mass, trim_with_profile,
};
use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
use crate::policy::NegLog;
use eyre::Report;
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
  pub fn resolve_posterior(
    &self,
    eps: f64,
    posterior_of: impl Fn(&Distribution<NegLog>) -> Result<Distribution<NegLog>, Report>,
  ) -> Result<ResolvedForward, Report> {
    let message = self
      .conv
      .resample(self.grid_points, self.max_points, OutputSpacing::Capped)?;
    let posterior = posterior_of(&message)?;
    let WindowedConvolution::Window {
      normalized,
      window: (lo, hi),
      ..
    } = &self.conv
    else {
      return Ok(ResolvedForward {
        message,
        posterior,
        spacing: MessageSpacing::Exact,
        profile: PosteriorProfile::Unknown,
      });
    };
    let resample = |max_dx: f64| {
      resample_to_mass_window(
        normalized,
        *lo,
        *hi,
        self.grid_points,
        self.max_points,
        SpacingLimit::AtMost(max_dx),
      )
    };

    let resolving_dx = match self.finer_spacing_needed(&message, &posterior, eps, normalized.dx()) {
      SpacingCheck::Resolved(profile) => {
        return Ok(ResolvedForward {
          message,
          posterior,
          spacing: MessageSpacing::Capped,
          profile,
        });
      },
      SpacingCheck::Finer(resolving_dx) => resolving_dx,
    };
    let message = resample(resolving_dx)?;
    let posterior = posterior_of(&message)?;
    if let SpacingCheck::Resolved(profile) = self.finer_spacing_needed(&message, &posterior, eps, normalized.dx()) {
      return Ok(ResolvedForward {
        message,
        posterior,
        spacing: MessageSpacing::Resolving,
        profile,
      });
    }

    let message = resample(normalized.dx())?;
    let posterior = posterior_of(&message)?;
    Ok(ResolvedForward {
      message,
      posterior,
      spacing: MessageSpacing::Fine,
      profile: PosteriorProfile::Unknown,
    })
  }

  fn finer_spacing_needed(
    &self,
    message: &Distribution<NegLog>,
    posterior: &Distribution<NegLog>,
    eps: f64,
    fine_dx: f64,
  ) -> SpacingCheck {
    let (Distribution::Function(message), Distribution::Function(posterior)) = (message, posterior) else {
      return SpacingCheck::Resolved(PosteriorProfile::Unknown);
    };
    let Some(profile) = mass_sizable_profile(posterior) else {
      return SpacingCheck::Resolved(PosteriorProfile::NotMassSizable);
    };
    let Some(fractions) = profile.fractions(TailMass::Required) else {
      return SpacingCheck::Resolved(PosteriorProfile::Unknown);
    };
    let mass_width = fractions.survival_inverse(eps) - fractions.quantile(eps);
    let intervals = self.grid_points.saturating_sub(1).max(1).to_f64_exact();
    let resolving_dx = fine_dx.max(mass_width / intervals);
    if message.dx() > resolving_dx {
      SpacingCheck::Finer(resolving_dx)
    } else {
      SpacingCheck::Resolved(PosteriorProfile::MassSizable(Box::new(profile)))
    }
  }
}

#[expect(
  clippy::partial_pub_fields,
  reason = "the mass profile of the posterior only serves its trimming inside this crate"
)]
pub struct ResolvedForward {
  pub message: Distribution<NegLog>,
  pub posterior: Distribution<NegLog>,
  pub spacing: MessageSpacing,
  profile: PosteriorProfile,
}

impl ResolvedForward {
  pub fn trim_posterior_to_mass(self, eps: f64) -> Result<Self, Report> {
    let Self {
      message,
      posterior,
      spacing,
      profile,
    } = self;
    let posterior = match profile {
      PosteriorProfile::Unknown => trim_to_mass(posterior, eps)?,
      PosteriorProfile::NotMassSizable => posterior,
      PosteriorProfile::MassSizable(profile) => trim_with_profile(posterior, &profile, eps)?,
    };
    Ok(Self {
      message,
      posterior,
      spacing,
      profile: PosteriorProfile::Unknown,
    })
  }
}

enum SpacingCheck {
  Finer(f64),
  Resolved(PosteriorProfile),
}

enum PosteriorProfile {
  Unknown,
  NotMassSizable,
  MassSizable(Box<MassProfile>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MessageSpacing {
  Exact,
  Capped,
  Resolving,
  Fine,
}

enum WindowedConvolution {
  Exact(Distribution<NegLog>),
  Window {
    normalized: DistributionFunction<f64, NegLog>,
    window: (f64, f64),
    cap: f64,
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
    let cap = coarser_operand_spacing(a, b);
    let is_fft = cap.is_some();
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
    let Some(cap) = cap else {
      return Ok(Self::Exact(Distribution::Function(conv)));
    };
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
      cap,
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
        let spacing = match spacing {
          OutputSpacing::Capped => SpacingLimit::AtMost(*cap),
          OutputSpacing::Fine => SpacingLimit::Input,
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
  let (Distribution::Function(a), Distribution::Function(b)) = (a, b) else {
    return None;
  };
  Some(a.dx().max(b.dx()))
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
