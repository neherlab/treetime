#[cfg(test)]
mod tests {
  //! Eigendecomposition invariants for GTR rate matrices.

  use crate::gtr::__tests__::generators::tests::generators::arb_gtr_nuc;
  use approx::assert_abs_diff_eq;
  use ndarray::Array2;
  use proptest::prelude::*;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// All eigenvalues of a valid rate matrix Q have non-positive real parts.
    /// This ensures the system is stable: probabilities don't explode as t -> inf.
    /// For reversible models, eigenvalues are real (Q is similar to a symmetric matrix).
    /// Positive eigenvalues would cause exp(lambda*t) -> inf, violating probability bounds.
    /// Formula: Re(eigvals[i]) <= 0 for all i
    #[test]
    fn test_prop_gtr_eigen_eigvals_nonpositive(gtr in arb_gtr_nuc()) {
      for (i, &eigval) in gtr.eigvals.iter().enumerate() {
        prop_assert!(
          eigval <= 1e-10,  // Allow small numerical error
          "Eigenvalue {i} = {eigval} is positive"
        );
      }
    }

    /// Exactly one eigenvalue equals zero, corresponding to the stationary
    /// distribution pi. The associated eigenvector is the all-ones vector (or pi
    /// for the left eigenvector). This zero eigenvalue ensures:
    ///
    ///   pi * Q = 0  (pi is left eigenvector with eigenvalue 0)
    ///   Q * 1 = 0   (1-vector is right eigenvector with eigenvalue 0)
    ///
    /// Multiple zero eigenvalues would indicate a reducible chain.
    #[test]
    fn test_prop_gtr_eigen_eigvals_one_zero_eigenvalue(gtr in arb_gtr_nuc()) {
      let zero_count = gtr.eigvals.iter().filter(|&&e| e.abs() < 1e-10).count();
      prop_assert!(
        zero_count == 1,
        "Expected exactly 1 zero eigenvalue, found {zero_count}: {:?}",
        gtr.eigvals
      );
    }

    /// The eigendecomposition Q = V * diag(lambda) * V^{-1} requires V * V^{-1} = I.
    /// This identity is essential for computing matrix exponentials:
    ///
    ///   exp(Qt) = V * diag(exp(lambda*t)) * V^{-1}
    ///
    /// Numerical errors in V or V^{-1} propagate to expQt, so this invariant
    /// validates the decomposition quality.
    #[test]
    fn test_prop_gtr_eigen_eigendecomposition_v_times_v_inv_is_identity(gtr in arb_gtr_nuc()) {
      let product = gtr.v.dot(&gtr.v_inv);
      assert_abs_diff_eq!(product, Array2::eye(4), epsilon = 1e-10);
    }
  }
}
