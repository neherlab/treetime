#[cfg(test)]
mod __tests__;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn jukes_cantor_distance(p: f64, n_states: usize) -> f64 {
  debug_assert!(
    n_states >= 2,
    "jukes_cantor_distance: n_states must be >= 2, got {n_states}"
  );

  let k = n_states as f64;
  let p_sat = (k - 1.0) / k;

  let p_max = p_sat * (1.0 - SATURATION_MARGIN);
  let p_clamped = p.clamp(0.0, p_max);

  if p_clamped == 0.0 {
    return 0.0;
  }

  -(k - 1.0) / k * (1.0 - k / (k - 1.0) * p_clamped).ln()
}

const SATURATION_MARGIN: f64 = 1e-6;
