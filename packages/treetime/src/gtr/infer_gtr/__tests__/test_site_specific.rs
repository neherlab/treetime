#[cfg(test)]
mod tests {
  use crate::gtr::__tests__::site_specific_support::simulate_counts;
  use crate::gtr::gtr_site_specific::{GTRSiteSpecific, GTRSiteSpecificParams};
  use crate::gtr::infer_gtr::site_specific::{
    InferGtrSiteSpecificOptions, build_gtr_site_specific, infer_gtr_site_specific_impl,
  };
  use approx::assert_abs_diff_eq;
  use ndarray::prelude::*;

  // TODO(investigate): circular inference validation
  //
  // These tests generate mutation counts (n_ija, T_ia, root_state) from the
  // model's own closed-form equation: n_ija[i,j,a] = pi[i,a] * W[i,j] * T[j,a] * mu[a].
  // The solver then recovers the planted parameters from these synthetic counts.
  // This validates solver self-consistency (convergence, normalization, update
  // equations) but not end-to-end correctness.
  //
  // What this catches:
  // - Solver convergence failures
  // - Wrong update equation signs or index transpositions
  // - Normalization bugs in W or pi updates
  // - Pseudocount regularization effects (bounded by pc/total_time)
  //
  // What this does NOT catch:
  // - Bugs in the mutation count extraction path (get_branch_mutation_matrix,
  //   accumulate_mutation_counts) which derives n_ija and T_ia from branch
  //   profiles during ancestral reconstruction
  // - Parent/child orientation mistakes in the counting path
  // - Interaction between site-specific GTR and sequence compression
  //
  // The golden-master test (test_gm_gtr_site_specific_infer) cross-validates
  // against v0's GTR_site_specific.infer() on the same synthetic counts at 1e-3.
  //
  // Missing: an end-to-end test that derives counts from a real tree+alignment
  // via the partition system, then compares inferred parameters against v0.
  // This requires partition integration (GTR -> enum/trait at partition call sites).

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

    // Check per-site pi recovery. With total_time=10000 and pc=1, pseudocount
    // distortion is ~0.01%, so 1e-3 tolerance is appropriate.
    for a in 0..2 {
      let inferred_pi = result.pi.column(a);
      let original_pi = gtr.pi.column(a);
      for i in 0..4 {
        assert_abs_diff_eq!(inferred_pi[i], original_pi[i], epsilon = 1e-3);
      }
    }

    // Check W recovery: relative rates should be approximately correct.
    // The absolute scale is absorbed into mu, so compare ratios.
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
    let p = inferred.expQt(0.5).unwrap();
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
