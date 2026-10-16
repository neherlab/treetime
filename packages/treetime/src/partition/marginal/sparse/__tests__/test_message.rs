#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::sparse::message::propagate_raw_per_site;
  use crate::partition::storage::sparse::{FixedProfiles, SparseSeqDistribution};
  use crate::partition::storage::var_pos_map::VarPosMap;
  use approx::assert_abs_diff_eq;
  use ndarray::{Array1, array, aview1};
  use treetime_primitives::{AsciiChar, LogLh};

  #[test]
  fn test_message_propagate_raw_per_site_forward() {
    let a = AsciiChar::from_byte_unchecked(b'A');
    let c = AsciiChar::from_byte_unchecked(b'C');
    let g = AsciiChar::from_byte_unchecked(b'G');

    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.site_rates = Some(Array1::from_vec(vec![
      0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.1,
    ]));

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

    let result = propagate_raw_per_site(&gtr, gtr.site_rates.as_ref().unwrap(), t, false, &seq_dis, None);

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

    assert_abs_diff_eq!(aview1(result.variable.get(0).unwrap().dis), expected_0, epsilon = 1e-14);
    assert_abs_diff_eq!(aview1(result.variable.get(5).unwrap().dis), expected_5, epsilon = 1e-14);
    assert_abs_diff_eq!(
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

    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.site_rates = Some(Array1::from_vec(vec![
      0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.1,
    ]));

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

    let result = propagate_raw_per_site(&gtr, gtr.site_rates.as_ref().unwrap(), t, true, &seq_dis, None);

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

    assert_abs_diff_eq!(aview1(result.variable.get(0).unwrap().dis), expected_0, epsilon = 1e-14);
    assert_abs_diff_eq!(aview1(result.variable.get(5).unwrap().dis), expected_5, epsilon = 1e-14);
    assert_abs_diff_eq!(
      aview1(result.variable.get(10).unwrap().dis),
      expected_10,
      epsilon = 1e-14
    );
  }
}
