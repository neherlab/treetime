#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::sparse::message::*;
  use crate::seq::composition::Composition;
  use approx::assert_abs_diff_eq;
  use ndarray::array;

  #[test]
  fn test_propagate_raw_per_site_forward() {
    let a = AsciiChar::from_byte_unchecked(b'A');
    let c = AsciiChar::from_byte_unchecked(b'C');
    let g = AsciiChar::from_byte_unchecked(b'G');

    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.set_site_rates(Array1::from_vec(vec![
      0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.1,
    ]));

    let t = 0.3;
    let variable = btreemap! {
      0_usize  => VarPos { dis: array![0.9, 0.03, 0.04, 0.03], state: a },
      5_usize  => VarPos { dis: array![0.1, 0.6,  0.2,  0.1 ], state: c },
      10_usize => VarPos { dis: array![0.05, 0.05, 0.8, 0.1 ], state: g },
    };

    let seq_dis = SparseSeqDistribution {
      variable,
      variable_indel: BTreeSet::new(),
      fixed: btreemap! {},
      fixed_counts: Composition::new(std::iter::empty::<AsciiChar>(), AsciiChar::from_byte_unchecked(b'-')),
      log_lh: LogLh::ZERO,
    };

    let result = propagate_raw_per_site(&gtr, t, false, &seq_dis, None);

    // JC69 closed-form (Jukes and Cantor 1969, DOI 10.1016/B978-1-4832-3211-9.50009-7):
    //   P(t)_ii = (1/4)(1 + 3*exp(-4/3 * mu * r * t))
    //   P(t)_ij = (1/4)(1 - exp(-4/3 * mu * r * t))
    // mu = 0.75 after GTR normalization (avg_transition = pi^T W pi), r = site rate.
    let expected_0 = array![
      0.8094601846762877,
      0.06064424518648728,
      0.06925132495073787,
      0.06064424518648728
    ];
    let expected_5 = array![
      0.16767825458589602,
      0.4420840726329092,
      0.22255941819529867,
      0.16767825458589602
    ];
    let expected_10 = array![
      0.05591089329029837,
      0.05591089329029837,
      0.7837450434516795,
      0.10443316996772378
    ];

    assert_abs_diff_eq!(result.variable[&0].dis, expected_0, epsilon = 1e-14);
    assert_abs_diff_eq!(result.variable[&5].dis, expected_5, epsilon = 1e-14);
    assert_abs_diff_eq!(result.variable[&10].dis, expected_10, epsilon = 1e-14);
  }

  #[test]
  fn test_propagate_raw_per_site_backward() {
    let a = AsciiChar::from_byte_unchecked(b'A');
    let c = AsciiChar::from_byte_unchecked(b'C');
    let g = AsciiChar::from_byte_unchecked(b'G');

    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.set_site_rates(Array1::from_vec(vec![
      0.5, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.1,
    ]));

    let t = 0.3;
    let variable = btreemap! {
      0_usize  => VarPos { dis: array![0.9, 0.03, 0.04, 0.03], state: a },
      5_usize  => VarPos { dis: array![0.1, 0.6,  0.2,  0.1 ], state: c },
      10_usize => VarPos { dis: array![0.05, 0.05, 0.8, 0.1 ], state: g },
    };

    let seq_dis = SparseSeqDistribution {
      variable,
      variable_indel: BTreeSet::new(),
      fixed: btreemap! {},
      fixed_counts: Composition::new(std::iter::empty::<AsciiChar>(), AsciiChar::from_byte_unchecked(b'-')),
      log_lh: LogLh::ZERO,
    };

    let result = propagate_raw_per_site(&gtr, t, true, &seq_dis, None);

    // JC69 P(t) is symmetric, so P(t)^T = P(t). Same expected values as forward.
    let expected_0 = array![
      0.8094601846762877,
      0.06064424518648728,
      0.06925132495073787,
      0.06064424518648728
    ];
    let expected_5 = array![
      0.16767825458589602,
      0.4420840726329092,
      0.22255941819529867,
      0.16767825458589602
    ];
    let expected_10 = array![
      0.05591089329029837,
      0.05591089329029837,
      0.7837450434516795,
      0.10443316996772378
    ];

    assert_abs_diff_eq!(result.variable[&0].dis, expected_0, epsilon = 1e-14);
    assert_abs_diff_eq!(result.variable[&5].dis, expected_5, epsilon = 1e-14);
    assert_abs_diff_eq!(result.variable[&10].dis, expected_10, epsilon = 1e-14);
  }
}
