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
  use crate::partition::optimize::dense::{PartitionContribution, get_coefficients};
  use approx::assert_abs_diff_eq;
  use ndarray::{Array1, Axis, array, concatenate};
  use rstest::rstest;

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

  fn test_gtr() -> GTR {
    jc69(JC69Params::default()).expect("JC69 creation failed")
  }

  fn dense_contribution(parent: Array1<f64>, child: Array1<f64>, gtr: &GTR) -> PartitionContribution {
    let parent_2d = parent.insert_axis(Axis(0));
    let child_2d = child.insert_axis(Axis(0));
    get_coefficients(&make_dense_seq_dis(parent_2d), &make_dense_seq_dis(child_2d), gtr)
  }

  fn sparse_site(
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

  #[rustfmt::skip]
  #[rstest]
  #[case::identical(       array![0.25, 0.25, 0.25, 0.25], array![0.25, 0.25, 0.25, 0.25])]
  #[case::skewed_parent(   array![0.9,  0.03, 0.03, 0.04], array![0.25, 0.25, 0.25, 0.25])]
  #[case::skewed_child(    array![0.25, 0.25, 0.25, 0.25], array![0.05, 0.05, 0.85, 0.05])]
  #[case::both_skewed(     array![0.7,  0.1,  0.1,  0.1 ], array![0.1,  0.1,  0.1,  0.7 ])]
  #[case::near_certain(    array![0.97, 0.01, 0.01, 0.01], array![0.01, 0.97, 0.01, 0.01])]
  #[case::one_dominant(    array![1.0,  0.0,  0.0,  0.0 ], array![0.25, 0.25, 0.25, 0.25])]
  #[trace]
  fn test_coefficient_invariant_nonneg_site_lh_at_zero(
    #[case] parent: Array1<f64>,
    #[case] child: Array1<f64>,
  ) {
    let gtr = test_gtr();
    let contribution = dense_contribution(parent, child, &gtr);
    let coeff_sum: f64 = contribution.coefficients.row(0).sum();
    assert!(coeff_sum >= 0.0, "Coefficient sum must be non-negative for valid probability messages, got {coeff_sum}");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::mult_2(  2.0)]
  #[case::mult_5(  5.0)]
  #[case::mult_10(10.0)]
  #[case::mult_100(100.0)]
  #[trace]
  fn test_coefficient_invariant_multiplicity_linearity(#[case] multiplicity: f64) {
    let gtr = test_gtr();
    let parent = array![0.6, 0.2, 0.1, 0.1];
    let child = array![0.3, 0.3, 0.2, 0.2];
    let branch_length = 0.05;

    let single = optimize::sparse::PartitionContribution {
      site_contributions: vec![sparse_site(&parent, &child, &gtr, 1.0)],
      gtr: gtr.clone(),
    };
    let single_metrics = evaluate_sparse_contribution(&single, branch_length).expect("valid branch length");

    let multi = optimize::sparse::PartitionContribution {
      site_contributions: vec![sparse_site(&parent, &child, &gtr, multiplicity)],
      gtr: gtr.clone(),
    };
    let multi_metrics = evaluate_sparse_contribution(&multi, branch_length).expect("valid branch length");

    assert_abs_diff_eq!(multi_metrics.log_lh.value(), multiplicity * single_metrics.log_lh.value(), epsilon = 1e-12);
    assert_abs_diff_eq!(multi_metrics.derivative, multiplicity * single_metrics.derivative, epsilon = 1e-12);
    assert_abs_diff_eq!(multi_metrics.second_derivative, multiplicity * single_metrics.second_derivative, epsilon = 1e-12);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rows_3(  3)]
  #[case::rows_10(10)]
  #[case::rows_50(50)]
  #[trace]
  fn test_coefficient_invariant_dense_sparse_equivalence(#[case] n_rows: usize) {
    let gtr = test_gtr();
    let parent = array![0.5, 0.2, 0.2, 0.1];
    let child = array![0.3, 0.4, 0.2, 0.1];
    let branch_length = 0.1;

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
      site_contributions: vec![sparse_site(&parent, &child, &gtr, n_rows as f64)],
      gtr: gtr.clone(),
    };
    let sparse_metrics = evaluate_sparse_contribution(&sparse_contrib, branch_length).expect("valid branch length");

    assert_abs_diff_eq!(dense_metrics.log_lh.value(), sparse_metrics.log_lh.value(), epsilon = 1e-10);
    assert_abs_diff_eq!(dense_metrics.derivative, sparse_metrics.derivative, epsilon = 1e-10);
    assert_abs_diff_eq!(dense_metrics.second_derivative, sparse_metrics.second_derivative, epsilon = 1e-10);
  }

  #[test]
  fn test_coefficient_invariant_additivity() {
    let gtr = test_gtr();
    let branch_length = 0.05;

    let parent_a = array![0.7, 0.1, 0.1, 0.1];
    let child_a = array![0.2, 0.6, 0.1, 0.1];
    let parent_b = array![0.1, 0.1, 0.7, 0.1];
    let child_b = array![0.1, 0.1, 0.1, 0.7];

    let contrib_a = dense_contribution(parent_a.clone(), child_a.clone(), &gtr);
    let contrib_b = dense_contribution(parent_b.clone(), child_b.clone(), &gtr);
    let metrics_a = evaluate_dense_contribution(&contrib_a, branch_length).expect("valid branch length");
    let metrics_b = evaluate_dense_contribution(&contrib_b, branch_length).expect("valid branch length");

    let parents = concatenate(
      Axis(0),
      &[
        parent_a.view().insert_axis(Axis(0)),
        parent_b.view().insert_axis(Axis(0)),
      ],
    )
    .unwrap();
    let children = concatenate(
      Axis(0),
      &[child_a.view().insert_axis(Axis(0)), child_b.view().insert_axis(Axis(0))],
    )
    .unwrap();
    let contrib_combined = get_coefficients(&make_dense_seq_dis(parents), &make_dense_seq_dis(children), &gtr);
    let metrics_combined = evaluate_dense_contribution(&contrib_combined, branch_length).expect("valid branch length");

    assert_abs_diff_eq!(
      metrics_combined.log_lh.value(),
      metrics_a.log_lh.value() + metrics_b.log_lh.value(),
      epsilon = 1e-12
    );
    assert_abs_diff_eq!(
      metrics_combined.derivative,
      metrics_a.derivative + metrics_b.derivative,
      epsilon = 1e-12
    );
    assert_abs_diff_eq!(
      metrics_combined.second_derivative,
      metrics_a.second_derivative + metrics_b.second_derivative,
      epsilon = 1e-12
    );
  }
}
