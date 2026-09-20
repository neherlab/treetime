#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::gtr::gtr::GTR;
  use crate::optimize::dense_eval::evaluate_dense_contribution;
  use crate::optimize::sparse_eval::evaluate_sparse_contribution;
  use crate::partition::optimize;
  use crate::partition::optimize::dense::get_coefficients;
  use ndarray::{Array1, Axis, array, concatenate};

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

  fn test_gtr() -> GTR {
    jc69(JC69Params::default()).expect("JC69 creation failed")
  }

  #[test]
  fn test_coefficient_boundary_disjoint_support_zero_at_t0() {
    let gtr = test_gtr();
    let parent = array![1.0, 0.0, 0.0, 0.0];
    let child = array![0.0, 1.0, 0.0, 0.0];

    let parent_2d = parent.insert_axis(Axis(0));
    let child_2d = child.insert_axis(Axis(0));
    let contribution = get_coefficients(&make_dense_seq_dis(parent_2d), &make_dense_seq_dis(child_2d), &gtr);

    let coeff_sum: f64 = contribution.coefficients.row(0).sum();
    assert!(
      coeff_sum.abs() < 1e-14,
      "Disjoint support should have ~zero coefficient sum at t=0, got {coeff_sum}"
    );

    let metrics = evaluate_dense_contribution(&contribution, 0.1).expect("valid branch length");
    assert!(
      metrics.log_lh.value().is_finite(),
      "At t>0, disjoint support should have finite log-lh"
    );
    assert!(metrics.derivative.is_finite(), "At t>0, derivative should be finite");
  }

  mod generators {
    use ndarray::Array1;
    use proptest::prelude::*;

    pub fn probability_vector() -> impl Strategy<Value = Array1<f64>> {
      prop::array::uniform4(1e-6..1.0_f64).prop_map(|raw| {
        let sum: f64 = raw.iter().sum();
        Array1::from_vec(raw.iter().map(|x| x / sum).collect())
      })
    }

    pub fn branch_length() -> impl Strategy<Value = f64> {
      1e-6..2.0_f64
    }

    pub fn multiplicity() -> impl Strategy<Value = f64> {
      1.0..500.0_f64
    }
  }

  mod prop_tests {
    use super::generators;
    use super::*;
    use crate::partition::optimize::dense::PartitionContribution;
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
    ) -> optimize::sparse::SiteContribution {
      let coefficients = parent.dot(&gtr.v) * child.dot(&gtr.v_inv.t());
      optimize::sparse::SiteContribution {
        multiplicity,
        coefficients,
      }
    }

    proptest! {
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

      #[test]
      fn test_prop_coefficient_multiplicity_linearity(
        parent in generators::probability_vector(),
        child in generators::probability_vector(),
        multiplicity in generators::multiplicity(),
        branch_length in generators::branch_length(),
      ) {
        let gtr = test_gtr();

        let single = optimize::sparse::PartitionContribution {
          site_contributions: vec![sparse_site_from(&parent, &child, &gtr, 1.0)],
          gtr: gtr.clone(),
        };
        let single_metrics = evaluate_sparse_contribution(&single, branch_length).expect("valid branch length");

        let multi = optimize::sparse::PartitionContribution {
          site_contributions: vec![sparse_site_from(&parent, &child, &gtr, multiplicity)],
          gtr: gtr.clone(),
        };
        let multi_metrics = evaluate_sparse_contribution(&multi, branch_length).expect("valid branch length");

        prop_assert_abs_diff_eq!(multi_metrics.log_lh.value(), multiplicity * single_metrics.log_lh.value(), epsilon = 1e-9);
        prop_assert_abs_diff_eq!(multi_metrics.derivative, multiplicity * single_metrics.derivative, epsilon = 1e-9);
        prop_assert_abs_diff_eq!(multi_metrics.second_derivative, multiplicity * single_metrics.second_derivative, epsilon = 1e-9);
      }

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
        let dense_metrics = evaluate_dense_contribution(&dense_contrib, branch_length).expect("valid branch length");

        let sparse_contrib = optimize::sparse::PartitionContribution {
          site_contributions: vec![sparse_site_from(&parent, &child, &gtr, n_rows as f64)],
          gtr: gtr.clone(),
        };
        let sparse_metrics = evaluate_sparse_contribution(&sparse_contrib, branch_length).expect("valid branch length");

        prop_assert_abs_diff_eq!(dense_metrics.log_lh.value(), sparse_metrics.log_lh.value(), epsilon = 1e-8);
        prop_assert_abs_diff_eq!(dense_metrics.derivative, sparse_metrics.derivative, epsilon = 1e-8);
        prop_assert_abs_diff_eq!(dense_metrics.second_derivative, sparse_metrics.second_derivative, epsilon = 1e-8);
      }

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
        let metrics_a = evaluate_dense_contribution(&contrib_a, branch_length).expect("valid branch length");
        let metrics_b = evaluate_dense_contribution(&contrib_b, branch_length).expect("valid branch length");

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
        let metrics_combined = evaluate_dense_contribution(&contrib_combined, branch_length).expect("valid branch length");

        prop_assert_abs_diff_eq!(metrics_combined.log_lh.value(), metrics_a.log_lh.value() + metrics_b.log_lh.value(), epsilon = 1e-9);
        prop_assert_abs_diff_eq!(metrics_combined.derivative, metrics_a.derivative + metrics_b.derivative, epsilon = 1e-9);
        prop_assert_abs_diff_eq!(metrics_combined.second_derivative, metrics_a.second_derivative + metrics_b.second_derivative, epsilon = 1e-9);
      }

      #[test]
      fn test_prop_coefficient_dense_finite_difference_derivative(
        parent in generators::probability_vector(),
        child in generators::probability_vector(),
        branch_length in 1e-4..1.0_f64,
      ) {
        let gtr = test_gtr();
        let contribution = dense_contribution_from(&parent, &child, &gtr);

        let metrics = evaluate_dense_contribution(&contribution, branch_length).expect("valid branch length");

        let h = f64::max(branch_length * 1e-4, 1e-5);
        let lh_plus = evaluate_dense_contribution(&contribution, branch_length + h).expect("valid branch length").log_lh.value();
        let lh_minus = evaluate_dense_contribution(&contribution, branch_length - h).expect("valid branch length").log_lh.value();
        let numerical_derivative = (lh_plus - lh_minus) / (2.0 * h);

        prop_assert_relative_eq!(metrics.derivative, numerical_derivative, max_relative = 1e-4);
      }

      #[test]
      fn test_prop_coefficient_dense_hessian_matches_d1_finite_difference(
        parent in generators::probability_vector(),
        child in generators::probability_vector(),
        branch_length in 1e-3..1.0_f64,
      ) {
        let gtr = test_gtr();
        let contribution = dense_contribution_from(&parent, &child, &gtr);

        let metrics = evaluate_dense_contribution(&contribution, branch_length).expect("valid branch length");

        let h = 1e-4;
        let d1_plus = evaluate_dense_contribution(&contribution, branch_length + h).expect("valid branch length").derivative;
        let d1_minus = evaluate_dense_contribution(&contribution, branch_length - h).expect("valid branch length").derivative;
        let numerical_second = (d1_plus - d1_minus) / (2.0 * h);

        prop_assert_relative_eq!(metrics.second_derivative, numerical_second, max_relative = 1e-5);
      }
    }
  }

  mod cancellation_regime {
    use super::*;
    use crate::optimize::sparse_eval::evaluate_sparse_contribution;
    use crate::partition::optimize::sparse::{PartitionContribution, SiteContribution};
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
      let site = SiteContribution {
        multiplicity: 1.0,
        coefficients: array![1.0, 0.0, 0.0, epsilon],
      };
      let contribution = PartitionContribution {
        site_contributions: vec![site],
        gtr,
      };

      let metrics = evaluate_sparse_contribution(&contribution, branch_length).expect("valid branch length");

      let lambda = -4.0 / 3.0;
      let exp_lt = (lambda * branch_length).exp();
      let s = exp_lt + epsilon;
      let w0 = exp_lt / s;
      let w3 = epsilon / s;
      let mean = w0 * lambda;
      let expected = w0 * (lambda - mean).powi(2) + w3 * mean * mean;

      let h = 1e-6;
      let d1_plus = evaluate_sparse_contribution(&contribution, branch_length + h).expect("valid branch length").derivative;
      let d1_minus = evaluate_sparse_contribution(&contribution, branch_length - h).expect("valid branch length").derivative;
      let numerical = (d1_plus - d1_minus) / (2.0 * h);

      assert_abs_diff_eq!(metrics.second_derivative, expected, epsilon = 1e-12);
      assert_abs_diff_eq!(metrics.second_derivative, numerical, epsilon = 1e-8);
    }
  }
}
