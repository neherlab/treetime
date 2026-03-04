#[cfg(test)]
mod tests {
  //! Transition matrix P(t) = exp(Qt) invariants.
  //!
  //! NOTE: This implementation uses a transposed convention where P(t) is
  //! column-stochastic (columns sum to 1), not row-stochastic.

  use crate::gtr::__tests__::generators::tests::generators::{arb_branch_len, arb_gtr_nuc, arb_profile_nuc};
  use crate::gtr::__tests__::prop_support::{prop_assert_columns_sum_to, prop_assert_rows_sum_to};
  use ndarray::Array2;
  use proptest::prelude::*;
  use treetime_utils::prop_assert_array_abs_diff_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// P(t) is a column-stochastic matrix in this transposed convention:
    /// each column is a probability distribution over destination states.
    ///
    ///   sum_i P[i,j](t) = 1 for all j
    ///
    /// This follows from probability conservation: Q has column sums = 0,
    /// so exp(Qt) has column sums = 1.
    ///
    /// Note: This is the TRANSPOSE of the standard textbook P(t) which has
    /// row sums = 1. The implementation uses P[i,j] = Pr(to state i | from j).
    #[test]
    fn test_prop_gtr_expqt_stochastic_columns(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      prop_assert_columns_sum_to(&p, 1.0, 1e-10)?;
    }

    /// All entries of P(t) are transition probabilities, hence non-negative:
    ///
    ///   P[i,j](t) = Pr(ending in state i | starting from state j) >= 0
    ///
    /// (Transposed convention: columns are source states, rows are destination states.)
    /// Negative "probabilities" would indicate numerical instability or an
    /// invalid rate matrix.
    #[test]
    fn test_prop_gtr_expqt_nonnegative(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      for i in 0..4 {
        for j in 0..4 {
          prop_assert!(
            p[[i, j]] >= -1e-14,  // Allow tiny negative due to clamping
            "P(t={t})[{i},{j}] = {} is negative",
            p[[i, j]]
          );
        }
      }
    }

    /// Transition probabilities cannot exceed 1:
    ///
    ///   P[i,j](t) <= 1 for all i, j, t
    ///
    /// Combined with non-negativity and row sums = 1, this fully constrains
    /// P(t) as a valid stochastic matrix.
    #[test]
    fn test_prop_gtr_expqt_bounded(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      for i in 0..4 {
        for j in 0..4 {
          prop_assert!(
            p[[i, j]] <= 1.0 + 1e-14,  // Allow tiny overshoot
            "P(t={t})[{i},{j}] = {} exceeds 1.0",
            p[[i, j]]
          );
        }
      }
    }

    /// At t = 0, no evolution has occurred, so the transition matrix is identity:
    ///
    ///   P(0) = exp(Q * 0) = exp(0) = I
    ///
    /// This means P[i,j](0) = 1 if i = j, else 0: the system starts in its
    /// initial state with probability 1.
    #[test]
    fn test_prop_gtr_expqt_zero_is_identity(gtr in arb_gtr_nuc()) {
      let p = gtr.expQt(0.0);
      prop_assert_array_abs_diff_eq!(p, Array2::eye(4), epsilon = 1e-10);
    }

    /// As t -> infinity, P(t) converges to the equilibrium distribution.
    /// In transposed convention where P[i,j] = Pr(to i | from j):
    ///
    ///   lim_{t -> inf} P[i,j](t) = pi[i]  for all j
    ///
    /// Each row i becomes constant at pi[i], regardless of starting state j.
    /// This is the ergodic theorem: the chain forgets its initial state.
    ///
    /// We scale time by 1000/mu to ensure sufficient evolution.
    ///
    /// Tolerance derivation (1e-6):
    /// - Effective time t_eff = 1000/mu with mu >= 0.001, so t_eff >= 100
    /// - GTR eigenvalues: one zero (stationary), others negative with |lambda| >= O(1)
    /// - Non-stationary contribution: exp(-|lambda| * t_eff) < exp(-100) ~ 1e-44
    /// - Remaining error from eigenvector reconstruction: V * diag(exp) * V_inv
    /// - Accumulated floating-point error in 4x4 matrix multiply: ~16 * eps ~ 4e-15
    /// - Conservative bound 1e-6 accounts for ill-conditioned eigenvector matrices
    #[test]
    fn test_prop_gtr_expqt_equilibrium_limit(gtr in arb_gtr_nuc()) {
      // Scale time by inverse of mu to ensure sufficient evolution
      let t = 1000.0 / gtr.mu.max(0.001);
      let p = gtr.expQt(t);
      // In transposed convention: P[i,j] -> pi[i] for all j
      // Each row of the equilibrium matrix is gtr.pi[i] repeated across columns
      let expected = Array2::from_shape_fn((4, 4), |(i, _j)| gtr.pi[i]);
      prop_assert_array_abs_diff_eq!(p, expected, epsilon = 1e-6);
    }

    /// The transition matrices form a semigroup under matrix multiplication:
    ///
    ///   P(s + t) = P(s) * P(t)
    ///
    /// This is the Chapman-Kolmogorov equation: the probability of transitioning
    /// over time s+t equals the sum over intermediate states of transitioning
    /// over s then t. Equivalently: exp(Q(s+t)) = exp(Qs) * exp(Qt).
    /// This property enables efficient likelihood computation on trees.
    #[test]
    fn test_prop_gtr_expqt_semigroup(gtr in arb_gtr_nuc(), s in 0.001_f64..1.0, t in 0.001_f64..1.0) {
      let p_s = gtr.expQt(s);
      let p_t = gtr.expQt(t);
      let p_st = gtr.expQt(s + t);
      let p_s_times_p_t = p_s.dot(&p_t);
      prop_assert_array_abs_diff_eq!(p_st, p_s_times_p_t, epsilon = 1e-10);
    }

    /// The equilibrium distribution pi is preserved under time evolution:
    /// In transposed convention, pi is the RIGHT eigenvector with eigenvalue 0:
    ///
    ///   P(t) @ pi = pi  for all t >= 0
    ///
    /// Starting from equilibrium, the system remains in equilibrium.
    /// This follows from Q @ pi = 0 (pi is right eigenvector of Q).
    #[test]
    fn test_prop_gtr_expqt_stationary_preserved(gtr in arb_gtr_nuc(), t in arb_branch_len()) {
      let p = gtr.expQt(t);
      let pi_evolved = p.dot(&gtr.pi);
      prop_assert_array_abs_diff_eq!(pi_evolved, gtr.pi, epsilon = 1e-10);
    }

    /// In ancestral reconstruction, two operations propagate information:
    ///
    ///   propagate_profile(p, t) = p @ P(t)      // child -> parent (backward)
    ///   evolve(p, t)            = p @ P(t).T    // parent -> child (forward)
    ///
    /// These differ by transpose because:
    /// - Backward: given child state distribution, what parent states explain it?
    /// - Forward: given parent state, what child states result from evolution?
    ///
    /// For a pure state e_i (one-hot vector), evolve gives column i of P(t),
    /// while propagate gives row i. The relationship:
    ///   evolve(p, t)[j] = sum_i p[i] * P[j,i](t)  (probability j generated p)
    ///   propagate(p, t)[j] = sum_i p[i] * P[i,j](t)  (probability p evolves to j)
    #[test]
    fn test_prop_gtr_expqt_evolve_transpose_of_propagate(
      gtr in arb_gtr_nuc(),
      profile in arb_profile_nuc(5),
      t in arb_branch_len()
    ) {
      let p = gtr.expQt(t);
      let p_t = p.t();

      // Manual computation for comparison
      let propagated = profile.dot(&p);
      let evolved = profile.dot(&p_t);

      // Compare with GTR methods
      let propagated_gtr = gtr.propagate_profile(&profile, t, false);
      let evolved_gtr = gtr.evolve(&profile, t, false);

      // Verify propagate_profile matches profile @ P(t)
      prop_assert_array_abs_diff_eq!(propagated_gtr, propagated, epsilon = 1e-10);

      // Verify evolve matches profile @ P(t).T
      prop_assert_array_abs_diff_eq!(evolved_gtr, evolved, epsilon = 1e-10);

      // In transposed convention, P is column-stochastic.
      // evolve uses P.T which is row-stochastic, so evolve preserves row sums.
      // propagate uses P which is column-stochastic, so propagate doesn't preserve row sums.
      prop_assert_rows_sum_to(&evolved_gtr, 1.0, 1e-10)?;
    }
  }
}
