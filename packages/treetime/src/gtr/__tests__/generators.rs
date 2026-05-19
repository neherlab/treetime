#[cfg(test)]
pub mod tests {
  pub mod generators {
    //! Proptest generators for GTR model parameters.
    //!
    //! Property-based testing (originating from Haskell's QuickCheck) defines logical
    //! properties that functions should satisfy, then generates random inputs to find
    //! counterexamples. Generators produce random values of specific types; when a test
    //! fails, the framework shrinks inputs to a minimal reproducing case.
    //!
    //! - <https://en.wikipedia.org/wiki/Software_testing#Property_testing>
    //! - <https://en.wikipedia.org/wiki/QuickCheck>
    //! - <https://docs.rs/proptest>

    use crate::gtr::gtr::{GTR, GTRParams};
    use ndarray::{Array1, Array2, Axis, stack};
    use proptest::prelude::*;

    /// Generate valid nucleotide equilibrium frequencies: positive, sum to 1.
    ///
    /// Uses staggered base values (0.5, 0.75, 1.0, 1.25) plus random offsets to ensure
    /// reasonable variation while avoiding extreme skew that causes numerical instability.
    /// This gives pi values in roughly [0.1, 0.4] range after normalization.
    pub fn arb_pi_nuc() -> impl Strategy<Value = Array1<f64>> {
      prop::collection::vec(0.0_f64..0.5, 4).prop_map(|offsets| {
        let bases = [0.5, 0.75, 1.0, 1.25];
        let raw: Vec<f64> = bases.iter().zip(offsets.iter()).map(|(b, o)| b + o).collect();
        let sum: f64 = raw.iter().sum();
        Array1::from_vec(raw.iter().map(|x| x / sum).collect())
      })
    }

    /// Generate valid amino acid equilibrium frequencies: positive, sum to 1.
    pub fn arb_pi_aa() -> impl Strategy<Value = Array1<f64>> {
      prop::collection::vec(0.01_f64..10.0, 20).prop_map(|raw| {
        let sum: f64 = raw.iter().sum();
        Array1::from_vec(raw.iter().map(|x| x / sum).collect())
      })
    }

    /// Generate valid nucleotide exchangeability matrix: symmetric, positive off-diagonal, zero diagonal.
    /// Upper triangle has 6 values for 4x4 matrix.
    ///
    /// Uses staggered base values (1.0, 1.5, 2.0, 2.5, 3.0, 3.5) plus random offsets to ensure
    /// eigenvalues are non-degenerate. Uniform W values create degenerate eigenspaces that
    /// cause numerical instability in the eigendecomposition.
    pub fn arb_w_nuc() -> impl Strategy<Value = Array2<f64>> {
      prop::collection::vec(0.0_f64..1.0, 6).prop_map(|offsets| {
        let bases = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5];
        let mut w = Array2::zeros((4, 4));
        let mut idx = 0;
        for i in 0..4 {
          for j in (i + 1)..4 {
            // Ensure non-uniform values: base + offset in [0, 1)
            w[[i, j]] = bases[idx] + offsets[idx];
            w[[j, i]] = w[[i, j]];
            idx += 1;
          }
        }
        w
      })
    }

    /// Generate valid amino acid exchangeability matrix: symmetric, positive off-diagonal, zero diagonal.
    /// Upper triangle has 190 values for 20x20 matrix.
    pub fn arb_w_aa() -> impl Strategy<Value = Array2<f64>> {
      prop::collection::vec(0.01_f64..10.0, 190).prop_map(|upper| {
        let mut w = Array2::zeros((20, 20));
        let mut idx = 0;
        for i in 0..20 {
          for j in (i + 1)..20 {
            w[[i, j]] = upper[idx];
            w[[j, i]] = upper[idx];
            idx += 1;
          }
        }
        w
      })
    }

    /// Generate branch lengths with log-uniform distribution.
    /// Covers both small (1e-10) and large (100) values.
    pub fn arb_branch_len() -> impl Strategy<Value = f64> {
      (-10.0_f64..2.0).prop_map(|exp| 10.0_f64.powf(exp))
    }

    /// Generate probability profile: L positions, each a probability distribution over 4 states.
    /// Uses range [0.1, 2.0] to avoid extreme skew that causes numerical instability.
    pub fn arb_profile_nuc(len: usize) -> impl Strategy<Value = Array2<f64>> {
      prop::collection::vec(prop::collection::vec(0.1_f64..2.0, 4), len).prop_map(|raw| {
        let rows: Vec<Array1<f64>> = raw
          .iter()
          .map(|row| {
            let sum: f64 = row.iter().sum();
            Array1::from_vec(row.iter().map(|x| x / sum).collect())
          })
          .collect();
        stack(Axis(0), &rows.iter().map(|r| r.view()).collect::<Vec<_>>()).expect("Stack should succeed")
      })
    }

    /// Generate a valid nucleotide GTR model.
    pub fn arb_gtr_nuc() -> impl Strategy<Value = GTR> {
      (arb_pi_nuc(), arb_w_nuc(), 0.1_f64..5.0).prop_map(|(pi, w, mu)| {
        GTR::new(GTRParams {
          n_states: 4,
          mu,
          W: Some(w),
          pi,
        })
        .expect("GTR construction should succeed with valid parameters")
      })
    }

    mod tests {
      //! Validate that generators produce outputs satisfying required invariants.
      //! Catches generator bugs at the source rather than as confusing numerical
      //! failures in downstream GTR tests.

      use super::*;
      use ndarray::{Array1, Axis};
      use treetime_utils::{
        prop_assert_abs_diff_eq, prop_assert_array_abs_diff_eq, prop_assert_array_diag_abs,
        prop_assert_array_offdiag_positive, prop_assert_array_positive,
      };

      proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        #[test]
        fn test_prop_generators_arb_pi_nuc_sums_to_one(pi in arb_pi_nuc()) {
          prop_assert_abs_diff_eq!(pi.sum(), 1.0, epsilon = 1e-14);
          prop_assert_array_positive!(pi);
        }

        #[test]
        fn test_prop_generators_arb_pi_aa_sums_to_one(pi in arb_pi_aa()) {
          prop_assert_abs_diff_eq!(pi.sum(), 1.0, epsilon = 1e-14);
          prop_assert_array_positive!(pi);
        }

        #[test]
        fn test_prop_generators_arb_w_nuc_symmetric_and_valid(w in arb_w_nuc()) {
          prop_assert_array_abs_diff_eq!(w, w.t().to_owned(), epsilon = 1e-14);
          prop_assert_array_diag_abs!(w, epsilon = 1e-14);
          prop_assert_array_offdiag_positive!(w);
        }

        #[test]
        fn test_prop_generators_arb_w_aa_symmetric_and_valid(w in arb_w_aa()) {
          prop_assert_array_abs_diff_eq!(w, w.t().to_owned(), epsilon = 1e-14);
          prop_assert_array_diag_abs!(w, epsilon = 1e-14);
        }

        #[test]
        fn test_prop_generators_arb_branch_len_range(t in arb_branch_len()) {
          prop_assert!((1e-10..=100.0).contains(&t), "branch length out of range: {t}");
          prop_assert!(t > 0.0, "branch length must be positive: {t}");
        }

        #[test]
        fn test_prop_generators_arb_profile_nuc_valid(profile in arb_profile_nuc(5)) {
          prop_assert_eq!(profile.shape(), &[5, 4]);
          let row_sums = profile.sum_axis(Axis(1));
          prop_assert_array_abs_diff_eq!(row_sums, Array1::ones(5), epsilon = 1e-14);
          prop_assert_array_positive!(profile);
        }

        #[test]
        fn test_prop_generators_arb_gtr_nuc_valid(gtr in arb_gtr_nuc()) {
          prop_assert_abs_diff_eq!(gtr.pi.sum(), 1.0, epsilon = 1e-10);
          prop_assert_array_abs_diff_eq!(gtr.W, gtr.W.t().to_owned(), epsilon = 1e-14);
          prop_assert!(gtr.mu > 0.0, "mu must be positive: {}", gtr.mu);
        }
      }
    }
  }
}
