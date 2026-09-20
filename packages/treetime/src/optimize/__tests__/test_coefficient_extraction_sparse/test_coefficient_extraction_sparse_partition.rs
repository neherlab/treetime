#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::optimize::sparse::{PartitionContribution, SiteContribution};
  use crate::pretty_assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_partition_contribution_empty_sites() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let contribution = PartitionContribution {
      site_contributions: vec![],
      gtr,
    };

    assert!(contribution.site_contributions.is_empty());
    assert_eq!(4, contribution.gtr.eigvals.len());
  }

  #[test]
  fn test_partition_contribution_single_variable_site() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let site = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.5, 0.2, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      gtr,
    };

    assert_eq!(1, contribution.site_contributions.len());
    pretty_assert_ulps_eq!(contribution.site_contributions[0].multiplicity, 1.0, max_ulps = 10);
  }

  #[test]
  fn test_partition_contribution_multiple_variable_sites() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let sites = vec![
      SiteContribution {
        multiplicity: 1.0,
        coefficients: array![0.5, 0.2, 0.2, 0.1],
      },
      SiteContribution {
        multiplicity: 1.0,
        coefficients: array![0.1, 0.4, 0.3, 0.2],
      },
      SiteContribution {
        multiplicity: 1.0,
        coefficients: array![0.3, 0.3, 0.2, 0.2],
      },
    ];

    let contribution = PartitionContribution {
      site_contributions: sites,
      gtr,
    };

    assert_eq!(3, contribution.site_contributions.len());
    for site in &contribution.site_contributions {
      pretty_assert_ulps_eq!(site.multiplicity, 1.0, max_ulps = 10);
    }
  }

  #[test]
  fn test_partition_contribution_fixed_sites_with_multiplicity() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let sites = vec![
      SiteContribution {
        multiplicity: 100.0,
        coefficients: array![0.9, 0.03, 0.03, 0.04],
      },
      SiteContribution {
        multiplicity: 50.0,
        coefficients: array![0.8, 0.1, 0.05, 0.05],
      },
    ];

    let contribution = PartitionContribution {
      site_contributions: sites,
      gtr,
    };

    assert_eq!(2, contribution.site_contributions.len());
    pretty_assert_ulps_eq!(contribution.site_contributions[0].multiplicity, 100.0, max_ulps = 10);
    pretty_assert_ulps_eq!(contribution.site_contributions[1].multiplicity, 50.0, max_ulps = 10);
  }
}
