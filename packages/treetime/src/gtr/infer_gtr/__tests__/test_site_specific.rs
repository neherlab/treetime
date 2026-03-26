#[cfg(test)]
mod tests {
  use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
  use crate::gtr::infer_gtr::site_specific::{
    InferGtrSiteSpecificOptions, MutationCountsSiteSpecific, build_gtr_site_specific, infer_gtr_site_specific_impl,
  };
  use approx::assert_abs_diff_eq;
  use ndarray::Array3;
  use ndarray::prelude::*;

  /// Simulate per-site mutation counts from a known site-specific GTR model.
  ///
  /// For each site a, generates synthetic n_ija and T_ia that are consistent
  /// with the model's W, pi_a, and mu_a. Uses the equation:
  ///   n_ija[i,j,a] = pi_a[i] * W[i,j] * T_ia[j,a] * mu_a
  fn simulate_counts(gtr: &GTRSiteSpecific, total_time: f64) -> MutationCountsSiteSpecific {
    let n = gtr.pi.nrows();
    let seq_len = gtr.seq_len;

    // T_ia: time in each state proportional to pi
    let mut T_ia = Array2::zeros((n, seq_len));
    for a in 0..seq_len {
      for i in 0..n {
        T_ia[[i, a]] = total_time * gtr.pi[[i, a]];
      }
    }

    // n_ija: expected mutation counts from the model
    let mut n_ija = Array3::zeros((n, n, seq_len));
    for a in 0..seq_len {
      for i in 0..n {
        for j in 0..n {
          if i != j {
            n_ija[[i, j, a]] = gtr.pi[[i, a]] * gtr.W[[i, j]] * T_ia[[j, a]] * gtr.mu[a];
          }
        }
      }
    }

    // root_state: proportional to pi
    let root_state = gtr.pi.clone();

    MutationCountsSiteSpecific {
      n_ija,
      T_ia,
      root_state,
    }
  }

  /// Inference from synthetic data should recover the original model's
  /// W matrix (shared) and per-site pi to reasonable accuracy.
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

    let counts = simulate_counts(&gtr, 1000.0);

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

    // Check per-site pi recovery (columns should match original pi proportions)
    for a in 0..2 {
      let inferred_pi = result.pi.column(a);
      let original_pi = gtr.pi.column(a);
      for i in 0..4 {
        assert_abs_diff_eq!(inferred_pi[i], original_pi[i], epsilon = 0.05);
      }
    }

    // Check W recovery: relative rates should be approximately correct
    // The absolute scale is absorbed into mu, so compare ratios
    let w_ref = result.W[[0, 1]];
    if w_ref > 1e-10 {
      let gtr_w_ref = gtr.W[[0, 1]];
      for i in 0..4 {
        for j in (i + 1)..4 {
          let inferred_ratio = result.W[[i, j]] / w_ref;
          let original_ratio = gtr.W[[i, j]] / gtr_w_ref;
          assert_abs_diff_eq!(inferred_ratio, original_ratio, epsilon = 0.15);
        }
      }
    }
  }

  /// The inferred model should produce valid transition matrices.
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

    // Verify the inferred model produces valid transition matrices
    let p = inferred.expQt(0.5);
    for a in 0..3 {
      let p_a = p.slice(s![.., .., a]);

      // Column-stochastic
      for j in 0..4 {
        let col_sum: f64 = p_a.column(j).sum();
        assert_abs_diff_eq!(col_sum, 1.0, epsilon = 1e-8);
      }

      // Non-negative
      for i in 0..4 {
        for j in 0..4 {
          assert!(
            p_a[[i, j]] >= -1e-14,
            "P[{i},{j},site={a}] = {} is negative",
            p_a[[i, j]]
          );
        }
      }
    }
  }
}
