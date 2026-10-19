use crate::policy::NegLog;
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use ndarray::Array1;
use treetime_grid::{BoundaryBehavior, Side, SoftTailLaw};
use treetime_utils::array::ndarray::min_or;
use treetime_utils::{make_error, make_internal_error};

pub fn distribution_multiply_by_fn<F, S>(
  distribution: &Distribution<NegLog>,
  weight: F,
  outward_slope: S,
) -> Result<Distribution<NegLog>, Report>
where
  F: Fn(f64) -> Result<f64, Report>,
  S: Fn(Side, f64) -> f64,
{
  if matches!(distribution, Distribution::Empty) {
    return Ok(Distribution::Empty);
  }
  if matches!(distribution, Distribution::Formula(_)) {
    return make_error!("distribution_multiply_by_fn requires a concrete Point, Range, or Function distribution");
  }

  let times = distribution.t();
  let (ordinates, amplitude) = match distribution {
    Distribution::Function(f) => (f.y().clone(), f.amplitude()),
    Distribution::Empty | Distribution::Point(..) | Distribution::Range(..) | Distribution::Formula(..) => {
      (distribution.y()?, 0.0)
    },
  };

  let weights: Array1<f64> = times.iter().map(|&t| weight(t)).collect::<Result<_, Report>>()?;
  let combined = ordinates + weights;

  let minimum = min_or(&combined, f64::INFINITY);
  if !minimum.is_finite() {
    return make_error!("distribution_multiply_by_fn found no finite weight over the distribution grid");
  }

  let result: Distribution<NegLog> = if let Distribution::Function(f) = distribution {
    let combined_fn = DistributionFunction::from_start_dx_values(f.x_min(), f.dx(), combined)?
      .with_left_extrap(weighted_tail_law(f, Side::Left, &outward_slope)?)?
      .with_right_extrap(weighted_tail_law(f, Side::Right, &outward_slope)?)?;
    Distribution::Function(combined_fn.with_amplitude(amplitude))
  } else {
    Distribution::function(times, combined)?
  };

  result.repeak()
}

fn weighted_tail_law(
  f: &DistributionFunction<f64, NegLog>,
  side: Side,
  outward_slope: impl Fn(Side, f64) -> f64,
) -> Result<BoundaryBehavior, Report> {
  let (behavior, edge) = match side {
    Side::Left => (f.left_extrap(), f.x_min()),
    Side::Right => (f.right_extrap(), f.x_max()),
  };
  let BoundaryBehavior::Linear(law) = behavior else {
    return Ok(behavior);
  };
  let weight_slope = outward_slope(side, edge);
  let composed = law.compose_multiply(&SoftTailLaw { slope: weight_slope });
  let decays = match side {
    Side::Left => composed.slope < 0.0,
    Side::Right => composed.slope > 0.0,
  };
  if !decays {
    return make_internal_error!(
      "The weight does not let the soft {side:?} tail decay: incoming slope {}, weight slope {weight_slope} beyond the edge at {edge}",
      law.slope
    );
  }
  Ok(BoundaryBehavior::Linear(composed))
}
