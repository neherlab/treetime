#[cfg(test)]
mod tests {
  use crate::partition::optimize_sparse::SiteContribution;
  use crate::pretty_assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_site_contribution_stores_multiplicity_and_coefficients() {
    let coefficients = array![0.1, 0.2, 0.3, 0.4];
    let contribution = SiteContribution {
      multiplicity: 5.0,
      coefficients: coefficients.clone(),
    };

    pretty_assert_ulps_eq!(contribution.multiplicity, 5.0, max_ulps = 10);
    for i in 0..4 {
      pretty_assert_ulps_eq!(contribution.coefficients[i], coefficients[i], max_ulps = 10);
    }
  }

  #[test]
  fn test_site_contribution_unit_multiplicity() {
    // Variable sites have multiplicity 1.0
    let contribution = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.25, 0.25, 0.25, 0.25],
    };

    pretty_assert_ulps_eq!(contribution.multiplicity, 1.0, max_ulps = 10);
  }
}
