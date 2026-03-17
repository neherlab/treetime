#[cfg(test)]
mod tests {
  //! Rate matrix Q invariants: structure, reversibility, and equilibrium frequencies.
  //!
  //! NOTE: This implementation uses a transposed Q convention where:
  //!   Q[i,j] = rate FROM state j TO state i (not the standard textbook convention)
  //!   Column sums = 0 (probability conservation per source state)
  //!   Detailed balance: pi[j] * Q[i,j] = pi[i] * Q[j,i]

  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::__tests__::generators::tests::generators::{arb_gtr_nuc, arb_pi_nuc, arb_w_nuc};
  use crate::gtr::__tests__::prop_support::{prop_assert_columns_sum_to, prop_assert_detailed_balance};
  use crate::gtr::gtr::{GTR, GTRParams};
  use proptest::prelude::*;
  use treetime_utils::{prop_assert_abs_diff_eq, prop_assert_array_abs_diff_eq};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// The rate matrix Q is the generator of a continuous-time Markov chain.
    /// In this transposed convention, each COLUMN must sum to zero because
    /// probability leaving source state j equals probability entering other states.
    ///
    /// This ensures probability conservation: d/dt(sum_i P[i,j](t)) = 0.
    /// Formula: sum(Q[:, j]) = 0 for all j
    #[test]
    fn test_prop_gtr_q_columns_sum_to_zero(gtr in arb_gtr_nuc()) {
      let q = gtr.Q();
      prop_assert_columns_sum_to(&q, 0.0, 1e-10)?;
    }

    /// Off-diagonal elements Q[i,j] (i != j) represent instantaneous rates
    /// of transitioning from state j to state i. Rates cannot be negative
    /// as they represent probability flow per unit time.
    /// Formula: Q[i,j] >= 0 for all i != j
    #[test]
    fn test_prop_gtr_q_offdiag_nonnegative(gtr in arb_gtr_nuc()) {
      let q = gtr.Q();
      for i in 0..4 {
        for j in 0..4 {
          if i != j {
            prop_assert!(
              q[[i, j]] >= 0.0,
              "Q[{i},{j}] = {} is negative",
              q[[i, j]]
            );
          }
        }
      }
    }

    /// Diagonal elements Q[j,j] represent the total rate of leaving state j.
    /// By convention, Q[j,j] = -sum_{i != j} Q[i,j], which is always <= 0.
    /// The magnitude |Q[j,j]| is the total exit rate from state j.
    /// Formula: Q[j,j] <= 0 for all j
    #[test]
    fn test_prop_gtr_q_diag_nonpositive(gtr in arb_gtr_nuc()) {
      let q = gtr.Q();
      for j in 0..4 {
        prop_assert!(
          q[[j, j]] <= 0.0,
          "Q[{j},{j}] = {} is positive",
          q[[j, j]]
        );
      }
    }

    /// A Markov chain is time-reversible if and only if detailed balance holds:
    /// the probability flux from j to i equals the flux from i to j at equilibrium.
    ///
    /// In transposed convention where Q[i,j] = rate from j to i:
    ///   pi[j] * Q[i,j] = pi[i] * Q[j,i]
    ///
    /// This ensures the process looks the same running forward or backward in time,
    /// which is required for valid phylogenetic likelihood calculations where we
    /// don't know the direction of evolution along branches.
    /// Reference: Kelly (1979), Reversibility and Stochastic Networks
    #[test]
    fn test_prop_gtr_q_detailed_balance(gtr in arb_gtr_nuc()) {
      let q = gtr.Q();
      prop_assert_detailed_balance(&q, &gtr.pi, 1e-10)?;
    }

    /// The exchangeability matrix W captures symmetric mutation rates independent
    /// of equilibrium frequencies. In this column-stochastic convention:
    ///
    ///   Q[i,j] = W[i,j] * pi[i]  (for i != j)
    ///
    /// W must be symmetric (W[i,j] = W[j,i]) to ensure detailed balance.
    /// This separates the intrinsic exchangeability between states from their
    /// equilibrium frequencies.
    #[test]
    fn test_prop_gtr_q_w_symmetric(gtr in arb_gtr_nuc()) {
      prop_assert_array_abs_diff_eq!(gtr.W, gtr.W.t().to_owned(), epsilon = 1e-14);
    }

    /// The equilibrium frequency vector pi represents the stationary distribution:
    /// the long-run proportion of time spent in each state. As a probability
    /// distribution, it must sum to 1.
    /// Formula: sum(pi) = 1
    #[test]
    fn test_prop_gtr_q_pi_sums_to_one(gtr in arb_gtr_nuc()) {
      prop_assert_abs_diff_eq!(gtr.pi.sum(), 1.0, epsilon = 1e-10);
    }

    /// Each equilibrium frequency must be strictly positive. A zero frequency
    /// would imply a state is never visited, making the rate matrix reducible.
    /// For irreducible chains (required for valid phylogenetic models), all
    /// states must be accessible.
    /// Formula: pi[i] > 0 for all i
    #[test]
    fn test_prop_gtr_q_pi_positive(gtr in arb_gtr_nuc()) {
      for (i, &p) in gtr.pi.iter().enumerate() {
        prop_assert!(
          p > 0.0,
          "pi[{i}] = {p} is not positive"
        );
      }
    }

    /// The overall rate mu scales time uniformly. A model with rate mu evolved
    /// for time t is equivalent to a model with rate 1 evolved for time mu*t:
    ///
    ///   exp(mu * Q_norm * t) = exp(Q_norm * (mu * t))
    ///
    /// This separates the substitution rate (mu) from the relative rates (Q_norm).
    /// In phylogenetics, branch lengths often absorb mu, measuring expected
    /// substitutions per site rather than time.
    ///
    /// Note: This test verifies that scaling mu by factor k and scaling t by 1/k
    /// produces the same transition matrix, within numerical tolerance.
    #[test]
    fn test_prop_gtr_q_mu_scaling((pi, w) in (arb_pi_nuc(), arb_w_nuc()), t in 0.01_f64..1.0) {
      let alphabet = Alphabet::new(AlphabetName::Nuc).expect("alphabet");
      let n_states = alphabet.n_canonical();

      // Create model with mu=1
      let gtr1 = GTR::new(GTRParams {
        n_states,
        mu: 1.0,
        W: Some(w.clone()),
        pi: pi.clone(),
      }).expect("GTR with mu=1");

      // Create model with mu=2
      let gtr2 = GTR::new(GTRParams {
        n_states,
        mu: 2.0,
        W: Some(w),
        pi,
      }).expect("GTR with mu=2");

      // P_1(2t) should equal P_2(t)
      let p1 = gtr1.expQt(2.0 * t);
      let p2 = gtr2.expQt(t);
      prop_assert_array_abs_diff_eq!(p1, p2, epsilon = 1e-10);
    }
  }
}
