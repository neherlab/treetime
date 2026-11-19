#[cfg(test)]
mod tests {
  use crate::gaussian::{
    GaussianParams, gaussian_convolution_pdf, gaussian_pdf, gaussian_pdf_grid, gaussian_product,
    gaussian_product_params,
  };
  use ndarray::{Array1, array};
  use std::f64::consts::PI;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[test]
  fn test_gaussian_product_params_of_one_gaussian_is_that_gaussian() {
    let result = gaussian_product_params(&[GaussianParams {
      mu: 1.0,
      sigma: 2.0,
      amplitude: 3.0,
    }]);

    pretty_assert_ulps_eq!(1.0, result.mu, max_ulps = 4);
    pretty_assert_ulps_eq!(2.0, result.sigma, max_ulps = 4);
    pretty_assert_ulps_eq!(3.0_f64.ln(), result.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_product_params_of_two_identical_gaussians_halves_the_variance() {
    let unit = GaussianParams {
      mu: 0.0,
      sigma: 1.0,
      amplitude: 1.0,
    };

    let result = gaussian_product_params(&[unit, unit]);

    pretty_assert_ulps_eq!(0.0, result.mu, max_ulps = 4);
    pretty_assert_ulps_eq!(1.0 / 2.0_f64.sqrt(), result.sigma, max_ulps = 4);
    pretty_assert_ulps_eq!(0.0, result.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_product_params_of_shifted_gaussians_scales_by_their_overlap() {
    let params = [
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

    pretty_assert_ulps_eq!(0.0, result.mu, max_ulps = 4);
    pretty_assert_ulps_eq!(1.0 / 2.0_f64.sqrt(), result.sigma, max_ulps = 4);
    pretty_assert_ulps_eq!(-1.0, result.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_product_params_of_unequal_widths_weights_means_by_precision() {
    let params = [
      GaussianParams {
        mu: 0.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
      GaussianParams {
        mu: 3.0,
        sigma: 2.0,
        amplitude: 1.0,
      },
    ];

    let result = gaussian_product_params(&params);

    pretty_assert_abs_diff_eq!(0.6, result.mu, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(0.8_f64.sqrt(), result.sigma, epsilon = 1e-15);
    pretty_assert_abs_diff_eq!(-0.5 * 9.0 / 5.0, result.log_scale, epsilon = 1e-15);
  }

  #[test]
  fn test_gaussian_product_params_of_no_gaussians_is_flat() {
    let result = gaussian_product_params(&[]);

    pretty_assert_ulps_eq!(0.0, result.mu, max_ulps = 4);
    assert!(result.sigma.is_infinite() && result.sigma > 0.0);
    pretty_assert_ulps_eq!(0.0, result.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_product_on_grid_is_the_unnormalized_product_gaussian() {
    let unit = GaussianParams {
      mu: 0.0,
      sigma: 1.0,
      amplitude: 1.0,
    };
    let grid = array![-2.0, -1.0, 0.0, 1.0, 2.0];

    let result = gaussian_product(&[unit, unit], &grid);

    let expected = grid.mapv(|x: f64| (-x * x).exp());
    pretty_assert_abs_diff_eq!(expected, result.normalized, epsilon = 1e-15);
    pretty_assert_ulps_eq!(0.0, result.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_product_of_no_gaussians_is_one_everywhere() {
    let result = gaussian_product(&[], &array![-1.0, 0.0, 2.0]);

    assert_eq!(Array1::<f64>::ones(3), result.normalized);
    pretty_assert_ulps_eq!(0.0, result.log_scale, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_pdf_matches_the_normal_density() {
    let expected = (-0.5 * 0.25_f64.powi(2)).exp() / (2.0 * (2.0 * PI).sqrt());

    pretty_assert_ulps_eq!(expected, gaussian_pdf(1.0, 2.0, 0.5), max_ulps = 4);
    pretty_assert_ulps_eq!(expected, gaussian_pdf(1.0, 2.0, 1.5), max_ulps = 4);
    pretty_assert_ulps_eq!(
      1.0 / (2.0 * (2.0 * PI).sqrt()),
      gaussian_pdf(1.0, 2.0, 1.0),
      max_ulps = 4
    );
  }

  #[test]
  fn test_gaussian_pdf_grid_matches_the_normal_density_at_each_point() {
    let result = gaussian_pdf_grid(0.0, 1.0, &array![-1.0, 0.0, 1.0]);

    let expected = array![(-0.5_f64).exp(), 1.0, (-0.5_f64).exp()] / (2.0 * PI).sqrt();
    pretty_assert_ulps_eq!(expected, result, max_ulps = 4);
  }

  #[test]
  fn test_gaussian_convolution_pdf_is_normal_with_summed_variances() {
    let peak = 1.0 / (2.0 * PI * 5.0).sqrt();

    pretty_assert_ulps_eq!(peak, gaussian_convolution_pdf(1.0, 2.0, 1.0, 1.0), max_ulps = 4);
    pretty_assert_ulps_eq!(
      peak * (-0.5 * 9.0 / 5.0_f64).exp(),
      gaussian_convolution_pdf(1.0, 2.0, 1.0, 4.0),
      max_ulps = 4
    );
  }
}
