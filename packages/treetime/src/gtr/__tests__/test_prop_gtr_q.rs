#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::__tests__::generators::tests::generators::{arb_gtr_nuc, arb_pi_nuc, arb_w_nuc};
  use crate::gtr::__tests__::prop_support::{prop_assert_columns_sum_to, prop_assert_detailed_balance};
  use crate::gtr::gtr::{GTR, GTRParams};
  use proptest::prelude::*;
  use treetime_utils::{
    prop_assert_abs_diff_eq, prop_assert_array_abs_diff_eq, prop_assert_array_diag_nonpositive,
    prop_assert_array_offdiag_nonneg, prop_assert_array_positive,
  };

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn test_prop_gtr_q_columns_sum_to_zero(gtr in arb_gtr_nuc()) {
      let q = gtr.Q();
      prop_assert_columns_sum_to(&q, 0.0, 1e-10)?;
    }

    #[test]
    fn test_prop_gtr_q_offdiag_nonnegative(gtr in arb_gtr_nuc()) {
      let q = gtr.Q();
      prop_assert_array_offdiag_nonneg!(q);
    }

    #[test]
    fn test_prop_gtr_q_diag_nonpositive(gtr in arb_gtr_nuc()) {
      let q = gtr.Q();
      prop_assert_array_diag_nonpositive!(q);
    }

    #[test]
    fn test_prop_gtr_q_detailed_balance(gtr in arb_gtr_nuc()) {
      let q = gtr.Q();
      prop_assert_detailed_balance(&q, &gtr.pi, 1e-10)?;
    }

    #[test]
    fn test_prop_gtr_q_w_symmetric(gtr in arb_gtr_nuc()) {
      prop_assert_array_abs_diff_eq!(gtr.W, gtr.W.t().to_owned(), epsilon = 1e-14);
    }

    #[test]
    fn test_prop_gtr_q_pi_sums_to_one(gtr in arb_gtr_nuc()) {
      prop_assert_abs_diff_eq!(gtr.pi.sum(), 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_prop_gtr_q_pi_positive(gtr in arb_gtr_nuc()) {
      prop_assert_array_positive!(gtr.pi);
    }

    #[test]
    fn test_prop_gtr_q_mu_scaling((pi, w) in (arb_pi_nuc(), arb_w_nuc()), t in 0.01_f64..1.0) {
      let alphabet = Alphabet::new(AlphabetName::Nuc).expect("alphabet");
      let n_states = alphabet.n_canonical();

      let gtr1 = GTR::new(GTRParams {
        n_states,
        mu: 1.0,
        W: Some(w.clone()),
        pi: pi.clone(),
      }).expect("GTR with mu=1");

      let gtr2 = GTR::new(GTRParams {
        n_states,
        mu: 2.0,
        W: Some(w),
        pi,
      }).expect("GTR with mu=2");

      let p1 = gtr1.expQt(2.0 * t);
      let p2 = gtr2.expQt(t);
      prop_assert_array_abs_diff_eq!(p1, p2, epsilon = 1e-10);
    }
  }
}
