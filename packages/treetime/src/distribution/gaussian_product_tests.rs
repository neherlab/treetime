#[cfg(test)]
mod tests {
  use crate::distribution::distribution::Distribution;
  use crate::distribution::distribution_function::DistributionFunction;
  use crate::distribution::scaled_distribution::ScaledDistribution;
  use crate::distribution::scaled_distribution_convolution::scaled_distribution_convolution;
  use crate::distribution::scaled_distribution_multiplication::{
    scaled_distribution_multiplication, scaled_distribution_multiply_many,
  };
  use crate::distribution::y_axis_policy::Plain;
  use approx::assert_relative_eq;
  use eyre::Report;
  use ndarray::Array1;

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
        assert_relative_eq!(f.y()[i], y[i], epsilon = 1e-10);
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
    assert_relative_eq!(product.peak_value(), 12.0, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gaussian_underflow_resistance_points() -> Result<(), Report> {
    let n = 100;
    let amplitude = 0.0001;

    let points: Vec<ScaledDistribution> =
      std::iter::repeat_with(|| create_point_scaled(0.0, amplitude)).take(n).collect();
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
    let g1 = create_gaussian_scaled(0.0, 1.0, 1.0, 201);
    let g2 = create_gaussian_scaled(0.0, 1.0, 1.0, 201);
    let product = scaled_distribution_multiplication(&g1, &g2)?;

    assert!(!product.is_empty());
    assert!(product.log_scale().is_finite());
    assert_relative_eq!(product.inner().max_value(), 1.0, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gaussian_product_ten() -> Result<(), Report> {
    let gaussians: Vec<ScaledDistribution> =
      std::iter::repeat_with(|| create_gaussian_scaled(0.0, 1.0, 1.0, 201)).take(10).collect();
    let refs: Vec<&ScaledDistribution> = gaussians.iter().collect();
    let product = scaled_distribution_multiply_many(&refs)?;

    assert!(!product.is_empty());
    assert!(product.log_scale().is_finite(), "log_scale must be finite for 10 Gaussians");
    Ok(())
  }

  #[test]
  fn test_gaussian_underflow_resistance_functions() -> Result<(), Report> {
    let n = 50;
    let amplitude = 0.01;

    let gaussians: Vec<ScaledDistribution> =
      std::iter::repeat_with(|| create_gaussian_scaled(0.0, 1.0, amplitude, 201)).take(n).collect();
    let refs: Vec<&ScaledDistribution> = gaussians.iter().collect();
    let product = scaled_distribution_multiply_many(&refs)?;

    assert!(!product.is_empty(), "Must not underflow to empty");
    assert!(product.log_scale().is_finite(), "Must handle small amplitudes without underflow");
    Ok(())
  }

  #[test]
  fn test_gaussian_convolution_basic() -> Result<(), Report> {
    let g1 = create_gaussian_scaled(0.0, 1.0, 1.0, 201);
    let g2 = create_gaussian_scaled(0.0, 1.0, 1.0, 201);

    let result = scaled_distribution_convolution(&g1, &g2)?;

    assert!(!result.is_empty());
    assert!(result.log_scale().is_finite());
    Ok(())
  }

  #[test]
  fn test_gaussian_analytical_formula() {
    struct GaussianParams {
      mu: f64,
      sigma: f64,
    }

    fn gaussian_product_analytical(params: &[GaussianParams]) -> (f64, f64, f64) {
      let precision_sum: f64 = params.iter().map(|p| 1.0 / p.sigma.powi(2)).sum();
      let sigma_star = (1.0 / precision_sum).sqrt();
      let mu_star = sigma_star.powi(2) * params.iter().map(|p| p.mu / p.sigma.powi(2)).sum::<f64>();

      let quadratic_term: f64 = params.iter().map(|p| (p.mu - mu_star).powi(2) / p.sigma.powi(2)).sum();
      let log_scale = -0.5 * quadratic_term;

      (mu_star, sigma_star, log_scale)
    }

    let params = vec![
      GaussianParams { mu: -1.0, sigma: 1.0 },
      GaussianParams { mu: 1.0, sigma: 1.0 },
    ];

    let (mu, sigma, log_scale) = gaussian_product_analytical(&params);

    assert_relative_eq!(mu, 0.0, epsilon = 1e-10);
    assert_relative_eq!(sigma, 1.0 / 2.0_f64.sqrt(), epsilon = 1e-10);
    assert_relative_eq!(log_scale, -1.0, epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_normalization_preserved() -> Result<(), Report> {
    let g1 = create_gaussian_scaled(0.0, 1.0, 10.0, 201);
    let g2 = create_gaussian_scaled(0.0, 1.0, 20.0, 201);
    let product = scaled_distribution_multiplication(&g1, &g2)?;

    assert!(!product.is_empty());
    assert_relative_eq!(product.inner().max_value(), 1.0, epsilon = 1e-10);
    Ok(())
  }
}
