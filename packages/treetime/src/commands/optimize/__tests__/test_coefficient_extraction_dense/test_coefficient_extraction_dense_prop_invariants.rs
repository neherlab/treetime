#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense::get_coefficients;
  use crate::commands::optimize::optimize_dense_eval::evaluate_dense_contribution;
  use crate::commands::optimize::optimize_sparse;
  use crate::commands::optimize::optimize_sparse_eval::evaluate_sparse_contribution;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use ndarray::{Array1, Axis, array, concatenate};

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

  fn test_gtr() -> GTR {
    jc69(JC69Params::default()).expect("JC69 creation failed")
  }

  // Boundary: disjoint support profiles have zero coefficient sum at t=0.
  // This is the singular case where ln(site_lh) = -inf at t=0.
  // The evaluator handles this mathematically (produces -inf log_lh)
  // but the optimizer must never reach this point (upstream guard in
  // run_optimize_mixed bumps BL away from zero).
  #[test]
  fn test_coefficient_boundary_disjoint_support_zero_at_t0() {
    let gtr = test_gtr();
    let parent = array![1.0, 0.0, 0.0, 0.0];
    let child = array![0.0, 1.0, 0.0, 0.0];

    let parent_2d = parent.insert_axis(Axis(0));
    let child_2d = child.insert_axis(Axis(0));
    let contribution = get_coefficients(&make_dense_seq_dis(parent_2d), &make_dense_seq_dis(child_2d), &gtr);

    // Coefficient sum at t=0 is zero for disjoint support
    let coeff_sum: f64 = contribution.coefficients.row(0).sum();
    assert!(
      coeff_sum.abs() < 1e-14,
      "Disjoint support should have ~zero coefficient sum at t=0, got {coeff_sum}"
    );

    // At positive branch length, site likelihood becomes positive
    // (transition matrix allows state changes)
    let metrics = evaluate_dense_contribution(&contribution, 0.1);
    assert!(
      metrics.log_lh.is_finite(),
      "At t>0, disjoint support should have finite log-lh"
    );
    assert!(metrics.derivative.is_finite(), "At t>0, derivative should be finite");
  }

  mod generators {
    use ndarray::Array1;
    use proptest::prelude::*;

    // Generate a valid probability vector of length 4 (nucleotide alphabet).
    // Uses the Dirichlet-like approach: draw 4 positive values, normalize.
    // Minimum component is 1e-6 to avoid exact zeros (which create the
    // disjoint-support singularity tested separately).
    pub fn probability_vector() -> impl Strategy<Value = Array1<f64>> {
      prop::array::uniform4(1e-6..1.0_f64).prop_map(|raw| {
        let sum: f64 = raw.iter().sum();
        Array1::from_vec(raw.iter().map(|x| x / sum).collect())
      })
    }

    // Generate a positive branch length spanning the biologically relevant range.
    pub fn branch_length() -> impl Strategy<Value = f64> {
      1e-6..2.0_f64
    }

    // Generate a multiplicity in a realistic range.
    pub fn multiplicity() -> impl Strategy<Value = f64> {
      1.0..500.0_f64
    }
  }

  mod prop_tests {
    use super::generators;
    use super::*;
    use crate::commands::optimize::optimize_dense::PartitionContribution;
    use proptest::prelude::*;
    use treetime_utils::{prop_assert_abs_diff_eq, prop_assert_relative_eq};

    fn dense_contribution_from(parent: &Array1<f64>, child: &Array1<f64>, gtr: &GTR) -> PartitionContribution {
      let parent_2d = parent.clone().insert_axis(Axis(0));
      let child_2d = child.clone().insert_axis(Axis(0));
      get_coefficients(&make_dense_seq_dis(parent_2d), &make_dense_seq_dis(child_2d), gtr)
    }

    fn sparse_site_from(
      parent: &Array1<f64>,
      child: &Array1<f64>,
      gtr: &GTR,
      multiplicity: f64,
    ) -> optimize_sparse::SiteContribution {
      let coefficients = parent.dot(&gtr.v) * child.dot(&gtr.v_inv.t());
      optimize_sparse::SiteContribution {
        multiplicity,
        coefficients,
      }
    }

    proptest! {
      // Invariant 1: non-negative site likelihood at t=0 for overlapping probability vectors.
      // With minimum component 1e-6, the inner product is always positive.
      #[test]
      fn test_prop_coefficient_nonneg_site_lh_at_zero(
        parent in generators::probability_vector(),
        child in generators::probability_vector(),
      ) {
        let gtr = test_gtr();
        let contribution = dense_contribution_from(&parent, &child, &gtr);
        let coeff_sum: f64 = contribution.coefficients.row(0).sum();
        prop_assert!(
          coeff_sum > 0.0,
          "Coefficient sum must be positive for overlapping probability vectors, got {coeff_sum}"
        );
      }

      // Invariant 2: multiplicity linearity.
      // sparse(multiplicity=m) == m * sparse(multiplicity=1) for all metrics.
      #[test]
      fn test_prop_coefficient_multiplicity_linearity(
        parent in generators::probability_vector(),
        child in generators::probability_vector(),
        multiplicity in generators::multiplicity(),
        branch_length in generators::branch_length(),
      ) {
        let gtr = test_gtr();

        let single = optimize_sparse::PartitionContribution {
          site_contributions: vec![sparse_site_from(&parent, &child, &gtr, 1.0)],
          gtr: gtr.clone(),
        };
        let single_metrics = evaluate_sparse_contribution(&single, branch_length);

        let multi = optimize_sparse::PartitionContribution {
          site_contributions: vec![sparse_site_from(&parent, &child, &gtr, multiplicity)],
          gtr: gtr.clone(),
        };
        let multi_metrics = evaluate_sparse_contribution(&multi, branch_length);

        prop_assert_abs_diff_eq!(multi_metrics.log_lh, multiplicity * single_metrics.log_lh, epsilon = 1e-9);
        prop_assert_abs_diff_eq!(multi_metrics.derivative, multiplicity * single_metrics.derivative, epsilon = 1e-9);
        prop_assert_abs_diff_eq!(multi_metrics.second_derivative, multiplicity * single_metrics.second_derivative, epsilon = 1e-9);
      }

      // Invariant 3: dense-sparse equivalence.
      // n identical dense rows == one sparse site with multiplicity n.
      #[test]
      fn test_prop_coefficient_dense_sparse_equivalence(
        parent in generators::probability_vector(),
        child in generators::probability_vector(),
        n_rows in 2..50_usize,
        branch_length in generators::branch_length(),
      ) {
        let gtr = test_gtr();

        let parent_2d = parent.view().insert_axis(Axis(0));
        let child_2d = child.view().insert_axis(Axis(0));
        let parents_stacked = concatenate(Axis(0), &vec![parent_2d; n_rows]).unwrap();
        let children_stacked = concatenate(Axis(0), &vec![child_2d; n_rows]).unwrap();
        let dense_contrib = get_coefficients(
          &make_dense_seq_dis(parents_stacked),
          &make_dense_seq_dis(children_stacked),
          &gtr,
        );
        let dense_metrics = evaluate_dense_contribution(&dense_contrib, branch_length);

        let sparse_contrib = optimize_sparse::PartitionContribution {
          site_contributions: vec![sparse_site_from(&parent, &child, &gtr, n_rows as f64)],
          gtr: gtr.clone(),
        };
        let sparse_metrics = evaluate_sparse_contribution(&sparse_contrib, branch_length);

        prop_assert_abs_diff_eq!(dense_metrics.log_lh, sparse_metrics.log_lh, epsilon = 1e-8);
        prop_assert_abs_diff_eq!(dense_metrics.derivative, sparse_metrics.derivative, epsilon = 1e-8);
        prop_assert_abs_diff_eq!(dense_metrics.second_derivative, sparse_metrics.second_derivative, epsilon = 1e-8);
      }

      // Invariant 4: coefficient additivity.
      // log_lh(site_a + site_b) == log_lh(site_a) + log_lh(site_b).
      #[test]
      fn test_prop_coefficient_additivity(
        parent_a in generators::probability_vector(),
        child_a in generators::probability_vector(),
        parent_b in generators::probability_vector(),
        child_b in generators::probability_vector(),
        branch_length in generators::branch_length(),
      ) {
        let gtr = test_gtr();

        let contrib_a = dense_contribution_from(&parent_a, &child_a, &gtr);
        let contrib_b = dense_contribution_from(&parent_b, &child_b, &gtr);
        let metrics_a = evaluate_dense_contribution(&contrib_a, branch_length);
        let metrics_b = evaluate_dense_contribution(&contrib_b, branch_length);

        let parents = concatenate(
          Axis(0),
          &[parent_a.view().insert_axis(Axis(0)), parent_b.view().insert_axis(Axis(0))],
        )
        .unwrap();
        let children = concatenate(
          Axis(0),
          &[child_a.view().insert_axis(Axis(0)), child_b.view().insert_axis(Axis(0))],
        )
        .unwrap();
        let contrib_combined = get_coefficients(&make_dense_seq_dis(parents), &make_dense_seq_dis(children), &gtr);
        let metrics_combined = evaluate_dense_contribution(&contrib_combined, branch_length);

        prop_assert_abs_diff_eq!(metrics_combined.log_lh, metrics_a.log_lh + metrics_b.log_lh, epsilon = 1e-9);
        prop_assert_abs_diff_eq!(metrics_combined.derivative, metrics_a.derivative + metrics_b.derivative, epsilon = 1e-9);
        prop_assert_abs_diff_eq!(metrics_combined.second_derivative, metrics_a.second_derivative + metrics_b.second_derivative, epsilon = 1e-9);
      }

      // Finite-difference derivative verification for dense evaluator.
      // Uses a fixed step with a floor to stay above the roundoff-dominated
      // regime: for `h < ~1e-8`, the `eps/h` rounding term dominates and the
      // central difference becomes unreliable as a derivative oracle.
      #[test]
      fn test_prop_coefficient_dense_finite_difference_derivative(
        parent in generators::probability_vector(),
        child in generators::probability_vector(),
        branch_length in 1e-4..1.0_f64,
      ) {
        let gtr = test_gtr();
        let contribution = dense_contribution_from(&parent, &child, &gtr);

        let metrics = evaluate_dense_contribution(&contribution, branch_length);

        let h = f64::max(branch_length * 1e-4, 1e-5);
        let lh_plus = evaluate_dense_contribution(&contribution, branch_length + h).log_lh;
        let lh_minus = evaluate_dense_contribution(&contribution, branch_length - h).log_lh;
        let numerical_derivative = (lh_plus - lh_minus) / (2.0 * h);

        prop_assert_relative_eq!(metrics.derivative, numerical_derivative, max_relative = 1e-4);
      }

      // Finite-difference verification of the analytical Hessian against the
      // central difference of the analytical FIRST derivative. This is far
      // more precise than the second difference of `log_lh` (which has
      // `O(eps / h^2)` rounding error in the numerator): a first-derivative
      // central difference has only `O(eps / h)` rounding error and `O(h^2)`
      // truncation error. With `h = 1e-4` the rounding floor is `~eps / h
      // ~ 1e-12` of the derivative magnitude, giving ~8-9 digits of
      // agreement -- well inside tolerance `max_relative = 1e-5`.
      //
      // Tight relative agreement here directly rules out catastrophic
      // cancellation in the analytical formula: the naive difference form
      // would miss it by 3-5 orders of magnitude on near-stationarity
      // inputs.
      #[test]
      fn test_prop_coefficient_dense_hessian_matches_d1_finite_difference(
        parent in generators::probability_vector(),
        child in generators::probability_vector(),
        branch_length in 1e-3..1.0_f64,
      ) {
        let gtr = test_gtr();
        let contribution = dense_contribution_from(&parent, &child, &gtr);

        let metrics = evaluate_dense_contribution(&contribution, branch_length);

        let h = 1e-4;
        let d1_plus = evaluate_dense_contribution(&contribution, branch_length + h).derivative;
        let d1_minus = evaluate_dense_contribution(&contribution, branch_length - h).derivative;
        let numerical_second = (d1_plus - d1_minus) / (2.0 * h);

        prop_assert_relative_eq!(metrics.second_derivative, numerical_second, max_relative = 1e-5);
      }
    }
  }

  // Regression guard for catastrophic cancellation in the Hessian formula.
  //
  // In the cancellation regime the posterior is concentrated at the nonzero
  // eigenvalue class, so both `E[lambda^2]` and `E[lambda]^2` approach the
  // same magnitude while their difference (the true variance) is small.
  // With the naive `E[lambda^2] - E[lambda]^2` form this regime loses up to
  // ~6 significant digits for `epsilon ~ 1e-6`; the centered (Welford) form
  // preserves full double-precision accuracy.
  //
  // The construction sets `k = [1, 0, 0, epsilon]` directly (bypassing
  // probability-vector coefficients so we can target the regime precisely).
  // `eigh` returns eigenvalues in ascending order, so for JC69 the eigvals
  // are `(-4/3, -4/3, -4/3, 0)` and this coefficient vector puts ~all
  // posterior mass on the `lambda = -4/3` class (index 0), making both
  // moments approach `16/9` while the variance is `O(epsilon)`.
  mod cancellation_regime {
    use super::*;
    use crate::commands::optimize::optimize_sparse::{PartitionContribution, SiteContribution};
    use crate::commands::optimize::optimize_sparse_eval::evaluate_sparse_contribution;
    use approx::assert_abs_diff_eq;
    use rstest::rstest;

    #[rustfmt::skip]
    #[rstest]
    #[case::eps_1e_2(1e-2)]
    #[case::eps_1e_3(1e-3)]
    #[case::eps_1e_4(1e-4)]
    #[case::eps_1e_5(1e-5)]
    #[case::eps_1e_6(1e-6)]
    #[case::eps_1e_8(1e-8)]
    #[trace]
    fn test_hessian_stable_in_cancellation_regime(#[case] epsilon: f64) {
      let gtr = test_gtr();
      let branch_length = 0.01;
      // k[0] sits on the lambda = -4/3 eigenvalue class (index 0 after
      // ascending sort); k[3] sits on the lambda = 0 class (index 3).
      let site = SiteContribution {
        multiplicity: 1.0,
        coefficients: array![1.0, 0.0, 0.0, epsilon],
      };
      let contribution = PartitionContribution {
        site_contributions: vec![site],
        gtr,
      };

      let metrics = evaluate_sparse_contribution(&contribution, branch_length);

      // Closed-form expected variance for this configuration.
      //
      //   `S      = e^{-4 t / 3} + epsilon`
      //   `w_0    = e^{-4 t / 3} / S`,  `w_3 = epsilon / S`,  `w_1 = w_2 = 0`
      //   `mean   = w_0 * (-4/3)`
      //   `var    = w_0 * (-4/3 - mean)^2 + w_3 * (0 - mean)^2`
      let lambda = -4.0 / 3.0;
      let exp_lt = (lambda * branch_length).exp();
      let s = exp_lt + epsilon;
      let w0 = exp_lt / s;
      let w3 = epsilon / s;
      let mean = w0 * lambda;
      let expected = w0 * (lambda - mean).powi(2) + w3 * mean * mean;

      // Cross-check the analytical Hessian against the central difference of
      // the analytical first derivative. The first-derivative central
      // difference has only `O(eps/h)` rounding error, so `h = 1e-6` gives
      // ~`1e-10` precision -- well below the assertion tolerance.
      let h = 1e-6;
      let d1_plus = evaluate_sparse_contribution(&contribution, branch_length + h).derivative;
      let d1_minus = evaluate_sparse_contribution(&contribution, branch_length - h).derivative;
      let numerical = (d1_plus - d1_minus) / (2.0 * h);

      assert_abs_diff_eq!(metrics.second_derivative, expected, epsilon = 1e-12);
      assert_abs_diff_eq!(metrics.second_derivative, numerical, epsilon = 1e-8);
    }
  }
}
