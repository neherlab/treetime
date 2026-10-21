#[cfg(test)]
pub(super) mod tests {
  use crate::partition::storage::dense::DenseSeqDistribution;
  use crate::partition::storage::sparse::SparseSeqDistribution;
  use crate::pretty_assert_ulps_eq;
  use ndarray::aview1;
  use treetime_utils::{pretty_assert_array_finite, pretty_assert_array_nonneg};

  pub(crate) fn assert_dense_profile_stable(profile: &DenseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.value().is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh.value()
    );
    pretty_assert_array_finite!(profile.dis);
    pretty_assert_array_nonneg!(profile.dis, epsilon = 1e-15);
    for row in profile.dis.outer_iter() {
      pretty_assert_ulps_eq!(row.sum(), 1.0, max_ulps = max_ulps);
    }
  }

  pub(crate) fn assert_sparse_profile_stable(profile: &SparseSeqDistribution, max_ulps: u32) {
    assert!(
      profile.log_lh.value().is_finite(),
      "Profile log_lh is not finite: {}",
      profile.log_lh.value()
    );

    for var_pos in profile.variable.values() {
      pretty_assert_array_finite!(aview1(var_pos.dis));
      pretty_assert_array_nonneg!(aview1(var_pos.dis), epsilon = 1e-15);
      pretty_assert_ulps_eq!(aview1(var_pos.dis).sum(), 1.0, max_ulps = max_ulps);
    }

    for fixed_dis in profile.fixed.iter().map(|fixed| &fixed.dis) {
      pretty_assert_array_finite!(aview1(fixed_dis));
      pretty_assert_array_nonneg!(aview1(fixed_dis), epsilon = 1e-15);
      pretty_assert_ulps_eq!(aview1(fixed_dis).sum(), 1.0, max_ulps = max_ulps);
    }
  }
}
