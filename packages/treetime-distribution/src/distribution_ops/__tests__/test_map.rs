#[cfg(test)]
mod tests {
  use crate::DistributionPlain;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_ops::map::distribution_map;
  use treetime_utils::assert_error;

  #[test]
  fn test_map_formula_returns_error() {
    let formula = DistributionPlain::Formula(DistributionFormula::new(|_| Ok(1.0), 0.0, 1.0));

    assert_error!(
      distribution_map(&formula, |value| value),
      "Cannot map Formula: operation not implemented"
    );
  }
}
