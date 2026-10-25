#[cfg(test)]
mod tests {
  use generators::{arb_gtr_with_profile, arb_time};
  use proptest::prelude::*;
  use treetime_utils::prop_assert_array_abs_diff_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn test_prop_gtr_propagate_profile_in_eigenbasis_matches_transition_matrix(
      (gtr, profile) in arb_gtr_with_profile(),
      t in arb_time(),
    ) {
      let expected = profile.dot(&gtr.expQt(t));
      let actual = gtr.propagate_profile(&profile, t, false);
      prop_assert_array_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    }

    #[test]
    fn test_prop_gtr_evolve_in_eigenbasis_matches_transition_matrix(
      (gtr, profile) in arb_gtr_with_profile(),
      t in arb_time(),
    ) {
      let expected = profile.dot(&gtr.expQt(t).t());
      let actual = gtr.evolve(&profile, t, false);
      prop_assert_array_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    }
  }

  mod generators {
    use crate::gtr::gtr::{GTR, Mu};
    use ndarray::Array2;
    use proptest::prelude::*;
    use std::cmp::Ordering;

    pub(super) fn arb_gtr_with_profile() -> impl Strategy<Value = (GTR, Array2<f64>)> {
      (2_usize..=120)
        .prop_flat_map(|n| (Just(n), prop_oneof![Just(1_usize), 1..n]))
        .prop_flat_map(|(n, rows)| {
          let pi = prop::collection::vec(0.01_f64..10.0, n);
          let exchangeabilities = prop::collection::vec(0.01_f64..10.0, n * n);
          let profile = prop::collection::vec(0.0_f64..1.0, rows * n);
          let mu = 0.1_f64..5.0;
          (pi, exchangeabilities, profile, mu).prop_map(move |(pi, exchangeabilities, profile, mu)| {
            let total: f64 = pi.iter().sum();
            let pi = pi.into_iter().map(|weight| weight / total).collect();
            let w = Array2::from_shape_fn((n, n), |(i, j)| match i.cmp(&j) {
              Ordering::Less => exchangeabilities[i * n + j],
              Ordering::Greater => exchangeabilities[j * n + i],
              Ordering::Equal => 0.0,
            });
            let gtr = GTR::builder()
              .n_states(n)
              .mu(Mu::FromMatrix(mu))
              .W(w)
              .pi(pi)
              .build()
              .expect("a GTR with positive frequencies and exchangeabilities is valid");
            let profile = Array2::from_shape_vec((rows, n), profile).expect("rows of n values");
            (gtr, profile)
          })
        })
    }

    pub(super) fn arb_time() -> impl Strategy<Value = f64> {
      (-6.0_f64..1.0).prop_map(|exponent| 10.0_f64.powf(exponent))
    }
  }
}
