#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_ops::map::distribution_map;
  use treetime_utils::assert_error;

  #[test]
  fn test_map_formula_returns_error() {
    // Oracle: kb/issues/H-distribution-result-api-panics-on-formula.md.
    let formula = Distribution::Formula(DistributionFormula::new(|_| Ok(1.0), 0.0, 1.0));

    assert_error!(
      distribution_map(&formula, |value| value),
      "Cannot map Formula: operation not implemented"
    );
  }
}
