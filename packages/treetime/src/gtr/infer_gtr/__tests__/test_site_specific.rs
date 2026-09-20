#[cfg(test)]
mod tests {
  use crate::gtr::__tests__::site_specific_support::simulate_counts;
  use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
  use crate::gtr::infer_gtr::site_specific::{
    InferGtrSiteSpecificOptions, build_gtr_site_specific, infer_gtr_site_specific_impl,
  };
  use approx::assert_abs_diff_eq;
  use ndarray::prelude::*;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_array_nonneg};

  #[test]
  fn test_infer_gtr_site_specific_recovers_parameters() {
    let pi = array![[0.1, 0.3], [0.2, 0.2], [0.3, 0.1], [0.4, 0.4]];
    let W = {
      let mut w = Array2::<f64>::ones([4, 4]);
      w[[0, 2]] = 3.0;
      w[[2, 0]] = 3.0;
      w[[1, 3]] = 2.0;
      w[[3, 1]] = 2.0;
      w
    };

    let gtr = GTRSiteSpecific::new(GTRSiteSpecificParams {
      n_states: 4,
      seq_len: 2,
      mu: array![1.0, 1.0],
      W: Some(W),
      pi,
      approximate: false,
    })
    .unwrap();

    let counts = simulate_counts(&gtr, 10000.0);

    let result = infer_gtr_site_specific_impl(
      &counts,
      &InferGtrSiteSpecificOptions {
        n_states: 4,
        pc: 1.0,
        max_iter: 50,
        dp: 1e-8,
        ..Default::default()
      },
    )
    .unwrap();

    for a in 0..2 {
      let inferred_pi = result.pi.column(a).to_owned();
      let original_pi = gtr.pi.column(a).to_owned();
      pretty_assert_abs_diff_eq!(inferred_pi, original_pi, epsilon = 1e-3);
    }

    let w_ref = result.W[[0, 1]];
    if w_ref > 1e-10 {
      let gtr_w_ref = gtr.W[[0, 1]];
      for i in 0..4 {
        for j in (i + 1)..4 {
          let inferred_ratio = result.W[[i, j]] / w_ref;
          let original_ratio = gtr.W[[i, j]] / gtr_w_ref;
          assert_abs_diff_eq!(inferred_ratio, original_ratio, epsilon = 1e-2);
        }
      }
    }
  }

  #[test]
  fn test_infer_gtr_site_specific_produces_valid_model() {
    let pi = array![[0.1, 0.25, 0.4], [0.2, 0.25, 0.1], [0.3, 0.25, 0.2], [0.4, 0.25, 0.3]];

    let gtr = GTRSiteSpecific::new(GTRSiteSpecificParams {
      n_states: 4,
      seq_len: 3,
      mu: array![1.0, 2.0, 0.5],
      W: None,
      pi,
      approximate: false,
    })
    .unwrap();

    let counts = simulate_counts(&gtr, 500.0);

    let result = infer_gtr_site_specific_impl(
      &counts,
      &InferGtrSiteSpecificOptions {
        n_states: 4,
        ..Default::default()
      },
    )
    .unwrap();

    let inferred = build_gtr_site_specific(&result, 4, false).unwrap();

    let p = inferred.expQt(0.5).unwrap();
    pretty_assert_array_nonneg!(p, epsilon = 1e-14);
    let col_sums = p.sum_axis(Axis(0));
    pretty_assert_abs_diff_eq!(col_sums, Array2::ones((4, 3)), epsilon = 1e-8);
  }
}
