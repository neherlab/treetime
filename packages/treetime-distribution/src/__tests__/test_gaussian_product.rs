#[cfg(test)]
mod tests {
  use crate::{
    Distribution, DistributionFunction, Plain, ScaledDistribution, scaled_distribution_convolution,
    scaled_distribution_multiplication, scaled_distribution_multiply_many,
  };
  use approx::{assert_relative_eq, assert_ulps_eq};
  use eyre::Report;
  use ndarray::Array1;
  use ordered_float::OrderedFloat;
  use treetime_analytical::{GaussianParams, gaussian_product_params};

  const GAUSSIAN_GRID_HALF_WIDTH_SIGMAS: f64 = 5.0;

  fn create_gaussian_scaled(mu: f64, sigma: f64, amplitude: f64, n_points: usize) -> ScaledDistribution {
    let x_min = mu - GAUSSIAN_GRID_HALF_WIDTH_SIGMAS * sigma;
    let x_max = mu + GAUSSIAN_GRID_HALF_WIDTH_SIGMAS * sigma;
    let dx = (x_max - x_min) / (n_points - 1) as f64;

    let y: Array1<f64> = Array1::from_shape_fn(n_points, |i| {
      let xi = x_min + dx * (i as f64);
      let exponent = -0.5 * ((xi - mu) / sigma).powi(2);
      amplitude * exponent.exp()
    });

    let dist_fn = DistributionFunction::from_start_dx_values(x_min, dx, y).unwrap();
    let dist = Distribution::<Plain>::Function(dist_fn);
    ScaledDistribution::from_plain(&dist)
  }

  fn create_point_scaled(t: f64, amplitude: f64) -> ScaledDistribution {
    let dist = Distribution::<Plain>::point(t, amplitude);
    ScaledDistribution::from_plain(&dist)
  }

  #[test]
  fn test_gaussian_scaled_distribution_round_trip() -> Result<(), Report> {
    let y = ndarray::array![0.1, 0.5, 1.0, 0.5, 0.1];
    let dist_fn = DistributionFunction::from_start_dx_values(0.0, 1.0, y.clone())?;
    let plain = Distribution::<Plain>::Function(dist_fn);

    let scaled = ScaledDistribution::from_plain(&plain);
    let recovered = scaled.to_plain();

    if let Distribution::Function(f) = recovered {
      for i in 0..5 {
        assert_ulps_eq!(f.y()[i], y[i], max_ulps = 4);
      }
    } else {
      panic!("Expected Function variant");
    }
    Ok(())
  }

  #[test]
  fn test_gaussian_point_multiplication_same_location() -> Result<(), Report> {
    let p1 = create_point_scaled(0.0, 3.0);
    let p2 = create_point_scaled(0.0, 4.0);
    let product = scaled_distribution_multiplication(&p1, &p2)?;

    assert!(!product.is_empty());
    assert_ulps_eq!(product.peak_value(), 12.0, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_gaussian_underflow_resistance_points() -> Result<(), Report> {
    let n = 100;
    let amplitude = 0.0001;

    let points: Vec<ScaledDistribution> = std::iter::repeat_with(|| create_point_scaled(0.0, amplitude))
      .take(n)
      .collect();
    let refs: Vec<&ScaledDistribution> = points.iter().collect();
    let product = scaled_distribution_multiply_many(&refs)?;

    assert!(!product.is_empty(), "Must not underflow to empty");
    assert!(product.log_scale().is_finite(), "Must handle 1e-400 without underflow");

    let expected_log_scale = (n as f64) * amplitude.ln();
    assert_relative_eq!(product.log_scale(), expected_log_scale, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gaussian_product_two_identical() -> Result<(), Report> {
    // Product of two N(0,1) Gaussians:
    // - sigma_star = 1/sqrt(2) ≈ 0.7071
    // - mu_star = 0
    // - log_scale = 0 (amplitudes=1, centered at mu_star)
    let g1 = create_gaussian_scaled(0.0, 1.0, 1.0, 201);
    let g2 = create_gaussian_scaled(0.0, 1.0, 1.0, 201);
    let product = scaled_distribution_multiplication(&g1, &g2)?;

    assert!(!product.is_empty());
    assert_ulps_eq!(product.inner().max_value(), 1.0, max_ulps = 4);

    // Verify log_scale is 0 (both amplitudes=1, both centered at result mean)
    assert_relative_eq!(product.log_scale(), 0.0, epsilon = 0.01);

    // Verify the width narrowed to sigma_star = 1/sqrt(2)
    // Width at half-max for Gaussian: FWHM = 2*sqrt(2*ln(2))*sigma ≈ 2.355*sigma
    // For sigma_star=0.7071, FWHM ≈ 1.665
    if let Distribution::Function(f) = product.inner() {
      let t = f.t();
      let y = f.y();
      let half_max_indices: Vec<usize> = y
        .iter()
        .enumerate()
        .filter(|&(_, v)| *v >= 0.49 && *v <= 0.51)
        .map(|(i, _)| i)
        .collect();
      if half_max_indices.len() >= 2 {
        let fwhm = t[*half_max_indices.last().unwrap()] - t[half_max_indices[0]];
        let expected_fwhm = 2.0 * (2.0 * 2.0_f64.ln()).sqrt() / 2.0_f64.sqrt();
        assert_relative_eq!(fwhm, expected_fwhm, epsilon = 0.1);
      }
    }
    Ok(())
  }

  #[test]
  fn test_gaussian_product_ten() -> Result<(), Report> {
    // Product of 10 N(0,1) Gaussians:
    // - sigma_star = 1/sqrt(10) ≈ 0.3162
    // - mu_star = 0
    // - log_scale = 0 (amplitudes=1, centered at mu_star)
    let gaussians: Vec<ScaledDistribution> = std::iter::repeat_with(|| create_gaussian_scaled(0.0, 1.0, 1.0, 201))
      .take(10)
      .collect();
    let refs: Vec<&ScaledDistribution> = gaussians.iter().collect();
    let product = scaled_distribution_multiply_many(&refs)?;

    assert!(!product.is_empty());

    // Verify log_scale is approximately 0
    assert_relative_eq!(product.log_scale(), 0.0, epsilon = 0.1);

    // Verify the width narrowed to sigma_star = 1/sqrt(10)
    if let Distribution::Function(f) = product.inner() {
      let t = f.t();
      let y = f.y();
      let half_max_indices: Vec<usize> = y
        .iter()
        .enumerate()
        .filter(|&(_, v)| *v >= 0.49 && *v <= 0.51)
        .map(|(i, _)| i)
        .collect();
      if half_max_indices.len() >= 2 {
        let fwhm = t[*half_max_indices.last().unwrap()] - t[half_max_indices[0]];
        let expected_fwhm = 2.0 * (2.0 * 2.0_f64.ln()).sqrt() / 10.0_f64.sqrt();
        assert_relative_eq!(fwhm, expected_fwhm, epsilon = 0.15);
      }
    }
    Ok(())
  }

  #[test]
  fn test_gaussian_underflow_resistance_functions() -> Result<(), Report> {
    let n = 50;
    let amplitude = 0.01;

    let gaussians: Vec<ScaledDistribution> =
      std::iter::repeat_with(|| create_gaussian_scaled(0.0, 1.0, amplitude, 201))
        .take(n)
        .collect();
    let refs: Vec<&ScaledDistribution> = gaussians.iter().collect();
    let product = scaled_distribution_multiply_many(&refs)?;

    assert!(!product.is_empty(), "Must not underflow to empty");
    assert!(
      product.log_scale().is_finite(),
      "Must handle small amplitudes without underflow"
    );
    Ok(())
  }

  #[test]
  fn test_gaussian_convolution_basic() -> Result<(), Report> {
    // Convolution of two N(0,1) Gaussians:
    // - sigma_conv = sqrt(1^2 + 1^2) = sqrt(2) ≈ 1.414
    // - mu_conv = 0
    let g1 = create_gaussian_scaled(0.0, 1.0, 1.0, 201);
    let g2 = create_gaussian_scaled(0.0, 1.0, 1.0, 201);

    let result = scaled_distribution_convolution(&g1, &g2)?;

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());

    // Verify the result is centered at 0 and has widened to sigma = sqrt(2)
    if let Distribution::Function(f) = result.inner() {
      let t = f.t();
      let y = f.y();

      // Find peak location (should be near 0)
      let max_idx = y
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .map(|(i, _)| i)
        .unwrap();
      let peak_t = t[max_idx];
      assert_relative_eq!(peak_t, 0.0, epsilon = 0.1);

      // Verify width increased to sigma = sqrt(2)
      // FWHM = 2*sqrt(2*ln(2))*sigma ≈ 2.355*sqrt(2) ≈ 3.33
      let half_max_indices: Vec<usize> = y
        .iter()
        .enumerate()
        .filter(|&(_, v)| *v >= 0.49 && *v <= 0.51)
        .map(|(i, _)| i)
        .collect();
      if half_max_indices.len() >= 2 {
        let fwhm = t[*half_max_indices.last().unwrap()] - t[half_max_indices[0]];
        let expected_fwhm = 2.0 * (2.0 * 2.0_f64.ln()).sqrt() * 2.0_f64.sqrt();
        assert_relative_eq!(fwhm, expected_fwhm, epsilon = 0.2);
      }
    }
    Ok(())
  }

  #[test]
  fn test_gaussian_analytical_formula() {
    let params = vec![
      GaussianParams {
        mu: -1.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
      GaussianParams {
        mu: 1.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
    ];

    let result = gaussian_product_params(&params);

    assert_ulps_eq!(result.mu, 0.0, max_ulps = 4);
    assert_ulps_eq!(result.sigma, 1.0 / 2.0_f64.sqrt(), max_ulps = 4);
    assert_ulps_eq!(result.log_scale, -1.0, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_product_matches_analytical() -> Result<(), Report> {
    // Two Gaussians at mu=-1 and mu=1, both sigma=1
    // Analytical: mu*=0, sigma*=1/sqrt(2), log_scale=-1.0
    let g1 = create_gaussian_scaled(-1.0, 1.0, 1.0, 401);
    let g2 = create_gaussian_scaled(1.0, 1.0, 1.0, 401);
    let product = scaled_distribution_multiplication(&g1, &g2)?;

    assert!(!product.is_empty());

    // Expected log_scale from analytical formula
    let expected_log_scale = -1.0;
    assert_relative_eq!(product.log_scale(), expected_log_scale, epsilon = 0.1);

    // Verify peak location is near mu*=0
    if let Distribution::Function(f) = product.inner() {
      let t = f.t();
      let y = f.y();
      let max_idx = y
        .iter()
        .enumerate()
        .max_by_key(|&(_, v)| OrderedFloat(*v))
        .map(|(i, _)| i)
        .unwrap();
      let peak_t = t[max_idx];
      assert_relative_eq!(peak_t, 0.0, epsilon = 0.05);
    }

    Ok(())
  }

  #[test]
  fn test_gaussian_normalization_preserved() -> Result<(), Report> {
    // Product of N(0,1) with amplitude=10 and N(0,1) with amplitude=20
    // Analytical: log_scale = ln(10) + ln(20) = ln(200) ≈ 5.298
    let g1 = create_gaussian_scaled(0.0, 1.0, 10.0, 201);
    let g2 = create_gaussian_scaled(0.0, 1.0, 20.0, 201);
    let product = scaled_distribution_multiplication(&g1, &g2)?;

    assert!(!product.is_empty());
    assert_ulps_eq!(product.inner().max_value(), 1.0, max_ulps = 4);

    // Verify log_scale captures the amplitude product
    let expected_log_scale = (10.0_f64 * 20.0).ln();
    assert_relative_eq!(product.log_scale(), expected_log_scale, epsilon = 0.1);

    // Verify peak value when unscaled: exp(log_scale) * max_value = 200
    let unscaled_peak = product.log_scale().exp() * product.inner().max_value();
    assert_relative_eq!(unscaled_peak, 200.0, epsilon = 1.0);
    Ok(())
  }
}
