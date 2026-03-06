use crate::commands::optimize::optimize_sparse::{PartitionContribution, SiteContribution};
use crate::gtr::get_gtr::{JC69Params, jc69};
use approx::assert_ulps_eq;
use ndarray::array;

#[test]
fn test_partition_contribution_empty_sites() {
  let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

  let contribution = PartitionContribution {
    site_contributions: vec![],
    eigenvalues: gtr.eigvals.to_owned(),
  };

  assert!(contribution.site_contributions.is_empty());
  assert_eq!(contribution.eigenvalues.len(), 4);
}

#[test]
fn test_partition_contribution_single_variable_site() {
  let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

  // One variable site with multiplicity 1
  let site = SiteContribution {
    multiplicity: 1.0,
    coefficients: array![0.5, 0.2, 0.2, 0.1],
  };

  let contribution = PartitionContribution {
    site_contributions: vec![site],
    eigenvalues: gtr.eigvals.to_owned(),
  };

  assert_eq!(contribution.site_contributions.len(), 1);
  assert_ulps_eq!(contribution.site_contributions[0].multiplicity, 1.0, max_ulps = 10);
}

#[test]
fn test_partition_contribution_multiple_variable_sites() {
  let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

  // Multiple variable sites, each with multiplicity 1
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
    eigenvalues: gtr.eigvals.to_owned(),
  };

  assert_eq!(contribution.site_contributions.len(), 3);
  for site in &contribution.site_contributions {
    assert_ulps_eq!(site.multiplicity, 1.0, max_ulps = 10);
  }
}

#[test]
fn test_partition_contribution_fixed_sites_with_multiplicity() {
  let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

  // Fixed sites have multiplicity > 1 (count of identical positions)
  let sites = vec![
    SiteContribution {
      multiplicity: 100.0, // 100 identical fixed positions
      coefficients: array![0.9, 0.03, 0.03, 0.04],
    },
    SiteContribution {
      multiplicity: 50.0, // 50 identical fixed positions
      coefficients: array![0.8, 0.1, 0.05, 0.05],
    },
  ];

  let contribution = PartitionContribution {
    site_contributions: sites,
    eigenvalues: gtr.eigvals.to_owned(),
  };

  assert_eq!(contribution.site_contributions.len(), 2);
  assert_ulps_eq!(contribution.site_contributions[0].multiplicity, 100.0, max_ulps = 10);
  assert_ulps_eq!(contribution.site_contributions[1].multiplicity, 50.0, max_ulps = 10);
}
