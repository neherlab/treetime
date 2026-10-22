pub fn within_epsilons_of_scale(a: f64, b: f64, scale: f64, epsilons: f64) -> bool {
  if a.is_finite() && b.is_finite() {
    (a - b).abs() <= epsilons * f64::EPSILON * scale
  } else {
    a.is_infinite() && b.is_infinite() && a.is_sign_positive() == b.is_sign_positive()
  }
}
