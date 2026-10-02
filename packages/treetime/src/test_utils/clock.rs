use treetime_utils::least_squares::LineFit;

pub(crate) fn half_residual_sum_of_squares(dates: &[f64], divs: &[f64]) -> f64 {
  let fit = LineFit::least_squares(dates, divs);
  0.5
    * dates
      .iter()
      .zip(divs)
      .map(|(date, div)| (div - fit.slope * date - fit.intercept).powi(2))
      .sum::<f64>()
}
