#[cfg(test)]
mod tests {
  use crate::partition::marginal::sparse::message::propagate_raw_per_site;
  use crate::partition::storage::sparse::{FixedProfiles, SparseSeqDistribution};
  use crate::partition::storage::var_pos_map::VarPosMap;
  use helpers::jc69_with_site_rates;
  use ndarray::{array, aview1};
  use treetime_primitives::{AsciiChar, LogLh};
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_message_propagate_raw_per_site_forward() {
    let a = AsciiChar::from_byte_unchecked(b'A');
    let c = AsciiChar::from_byte_unchecked(b'C');
    let g = AsciiChar::from_byte_unchecked(b'G');

    let gtr = jc69_with_site_rates(array![0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.1,]);

    let t = 0.3;
    #[rustfmt::skip]
    let sites = [
      (0,  a, [0.9,  0.03, 0.04, 0.03]),
      (5,  c, [0.1,  0.6,  0.2,  0.1 ]),
      (10, g, [0.05, 0.05, 0.8,  0.1 ]),
    ];
    let mut variable = VarPosMap::default();
    for (pos, state, dis) in sites {
      variable.insert(pos, state, &dis);
    }

    let seq_dis = SparseSeqDistribution {
      variable,
      fixed: FixedProfiles::default(),
      log_lh: LogLh::ZERO,
    };

    let result = propagate_raw_per_site(&gtr, gtr.site_rates().unwrap(), t, false, &seq_dis);

    let expected_0 = array![
      0.7821749895006883,
      0.06987923432284399,
      0.07806654185362383,
      0.06987923432284399
    ];
    let expected_5 = array![
      0.18260065538241677,
      0.40726513744102766,
      0.22753355179413895,
      0.18260065538241677
    ];
    let expected_10 = array![
      0.05784211216953536,
      0.05784211216953536,
      0.7784341915337779,
      0.10588158412715153
    ];

    pretty_assert_abs_diff_eq!(aview1(result.variable.get(0).unwrap().dis), expected_0, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(aview1(result.variable.get(5).unwrap().dis), expected_5, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(
      aview1(result.variable.get(10).unwrap().dis),
      expected_10,
      epsilon = 1e-14
    );
  }

  #[test]
  fn test_message_propagate_raw_per_site_backward() {
    let a = AsciiChar::from_byte_unchecked(b'A');
    let c = AsciiChar::from_byte_unchecked(b'C');
    let g = AsciiChar::from_byte_unchecked(b'G');

    let gtr = jc69_with_site_rates(array![0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.1,]);

    let t = 0.3;
    #[rustfmt::skip]
    let sites = [
      (0,  a, [0.9,  0.03, 0.04, 0.03]),
      (5,  c, [0.1,  0.6,  0.2,  0.1 ]),
      (10, g, [0.05, 0.05, 0.8,  0.1 ]),
    ];
    let mut variable = VarPosMap::default();
    for (pos, state, dis) in sites {
      variable.insert(pos, state, &dis);
    }

    let seq_dis = SparseSeqDistribution {
      variable,
      fixed: FixedProfiles::default(),
      log_lh: LogLh::ZERO,
    };

    let result = propagate_raw_per_site(&gtr, gtr.site_rates().unwrap(), t, true, &seq_dis);

    let expected_0 = array![
      0.7821749895006883,
      0.06987923432284399,
      0.07806654185362383,
      0.06987923432284399
    ];
    let expected_5 = array![
      0.18260065538241677,
      0.40726513744102766,
      0.22753355179413895,
      0.18260065538241677
    ];
    let expected_10 = array![
      0.05784211216953536,
      0.05784211216953536,
      0.7784341915337779,
      0.10588158412715153
    ];

    pretty_assert_abs_diff_eq!(aview1(result.variable.get(0).unwrap().dis), expected_0, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(aview1(result.variable.get(5).unwrap().dis), expected_5, epsilon = 1e-14);
    pretty_assert_abs_diff_eq!(
      aview1(result.variable.get(10).unwrap().dis),
      expected_10,
      epsilon = 1e-14
    );
  }

  mod helpers {
    use crate::gtr::gtr::{GTR, Mu};
    use ndarray::Array1;

    pub(super) fn jc69_with_site_rates(site_rates: Array1<f64>) -> GTR {
      GTR::builder()
        .n_states(4)
        .mu(Mu::Normalized(1.0))
        .pi(Array1::ones(4))
        .site_rates(site_rates)
        .build()
        .unwrap()
    }
  }
}
