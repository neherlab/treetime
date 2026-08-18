/// Threshold below which the ordinary-least-squares design matrix is treated as singular.
///
/// `denom = n*sum(x^2) - (sum x)^2` is `n^2` times the variance of the abscissae, so it vanishes
/// when every `x` is equal. Below this magnitude the slope is numerically undefined and the fit
/// falls back to a flat line through the mean ordinate.
const SINGULAR_DENOM_EPS: f64 = 1e-30;

/// A simple least-squares line fit `y = slope * x + intercept`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineFit {
  pub slope: f64,
  pub intercept: f64,
}

impl LineFit {
  /// Fit the least-squares line `y = slope * x + intercept` to paired samples.
  ///
  /// Solves the ordinary-least-squares normal equations for the two parameters. When the abscissae
  /// are degenerate -- every `x` equal, so the design matrix is singular -- the slope is undefined
  /// and the fit falls back to a flat line through the mean ordinate (`slope = 0`,
  /// `intercept = mean(y)`).
  ///
  /// `xs` and `ys` must have the same non-zero length.
  pub fn least_squares(xs: &[f64], ys: &[f64]) -> LineFit {
    let n = xs.len() as f64;
    let sum_x: f64 = xs.iter().sum();
    let sum_y: f64 = ys.iter().sum();
    let sum_xx: f64 = xs.iter().map(|x| x * x).sum();
    let sum_xy: f64 = xs.iter().zip(ys).map(|(x, y)| x * y).sum();

    let sum_x_squared = sum_x * sum_x;
    let denom = n * sum_xx - sum_x_squared;
    if denom.abs() < SINGULAR_DENOM_EPS {
      return LineFit {
        slope: 0.0,
        intercept: sum_y / n,
      };
    }

    let slope = (n * sum_xy - sum_x * sum_y) / denom;
    let intercept = (sum_y - slope * sum_x) / n;
    LineFit { slope, intercept }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use rstest::rstest;

  // Oracle: each case is an exact line `y = slope*x + intercept` sampled at the given abscissae, so
  // the fit must recover the generating parameters. The degenerate case has a singular design
  // matrix (all `x` equal) and falls back to `slope = 0`, `intercept = mean(y)`.
  #[rustfmt::skip]
  #[rstest]
  #[case::unit_slope(   &[0.0, 1.0, 2.0, 3.0], &[1.0, 3.0,  5.0,  7.0], ( 2.0, 1.0))]
  #[case::negative(     &[0.0, 1.0, 2.0, 3.0], &[5.0, 3.0,  1.0, -1.0], (-2.0, 5.0))]
  #[case::flat(         &[0.0, 1.0, 2.0, 3.0], &[4.0, 4.0,  4.0,  4.0], ( 0.0, 4.0))]
  #[case::offset_x(     &[10.0, 11.0, 12.0],   &[0.0, 0.5,  1.0],       ( 0.5, -5.0))]
  #[case::degenerate_x( &[2.0, 2.0, 2.0],      &[1.0, 2.0,  3.0],       ( 0.0, 2.0))]
  #[trace]
  fn test_least_squares_fit(#[case] xs: &[f64], #[case] ys: &[f64], #[case] expected: (f64, f64)) {
    let (slope, intercept) = expected;
    let fit = LineFit::least_squares(xs, ys);
    assert_abs_diff_eq!(slope, fit.slope, epsilon = 1e-12);
    assert_abs_diff_eq!(intercept, fit.intercept, epsilon = 1e-12);
  }
}
