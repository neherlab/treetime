#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::helpers::three_tip_tree;
  use crate::branch_lengths::{exact, perturbed, realized};
  use crate::substitution_model::SubstitutionModel;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::Seq;
  use treetime_utils::assert_error;
  use treetime_utils::sync::random::get_random_number_generator;

  #[test]
  fn test_branch_lengths_exact_is_rate_times_duration() {
    let tree = three_tip_tree();
    let model = SubstitutionModel::jc69(0.5).unwrap();
    let durations = btreemap! { tree.root_inner => 0.25, tree.root_c => 1.5, tree.inner_a => 0.0, tree.inner_b => 3.0 };
    let expected =
      btreemap! { tree.root_inner => 0.125, tree.root_c => 0.75, tree.inner_a => 0.0, tree.inner_b => 1.5 };
    let actual = exact(&durations, &model).unwrap();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(-1.0,          "Duration of edge 1 must be finite and non-negative, got -1")]
  #[case::nan(     f64::NAN,      "Duration of edge 1 must be finite and non-negative, got NaN")]
  #[case::infinite(f64::INFINITY, "Duration of edge 1 must be finite and non-negative, got inf")]
  #[trace]
  fn test_branch_lengths_exact_rejects_invalid_durations(#[case] duration: f64, #[case] expected: &str) {
    let tree = three_tip_tree();
    let model = SubstitutionModel::jc69(0.5).unwrap();
    let durations = btreemap! { tree.root_inner => 0.25, tree.root_c => duration, tree.inner_a => 0.0, tree.inner_b => 3.0 };
    let result = exact(&durations, &model);
    assert_error!(result, expected);
  }

  #[test]
  fn test_branch_lengths_exact_rejects_a_length_that_overflows() {
    let tree = three_tip_tree();
    let model = SubstitutionModel::jc69(2.0).unwrap();
    let durations =
      btreemap! { tree.root_inner => 0.25, tree.root_c => f64::MAX, tree.inner_a => 0.0, tree.inner_b => 3.0 };
    let result = exact(&durations, &model);
    assert_error!(
      result,
      "Branch length of edge 1 is not finite: substitution rate (mu) 2e0 times duration 1.7976931348623157e308"
    );
  }

  #[test]
  fn test_branch_lengths_realized_counts_differing_sites_per_site() {
    let tree = three_tip_tree();
    let sequences = btreemap! {
      tree.root  => Seq::try_from_str("ACGT").unwrap(),
      tree.inner => Seq::try_from_str("ACGA").unwrap(),
      tree.tip_c => Seq::try_from_str("TCGT").unwrap(),
      tree.tip_a => Seq::try_from_str("ACGA").unwrap(),
      tree.tip_b => Seq::try_from_str("GCTA").unwrap(),
    };
    let expected = btreemap! { tree.root_inner => 0.25, tree.root_c => 0.25, tree.inner_a => 0.0, tree.inner_b => 0.5 };
    let actual = realized(&tree.graph, &sequences).unwrap();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::different_lengths("ACGT", "ACG", "Sequences of nodes 0 and 1 have different lengths: 4 and 3")]
  #[case::empty(            "",     "",    "Sequence of node 0 is empty")]
  #[trace]
  fn test_branch_lengths_realized_rejects_invalid_sequences(
    #[case] root: &str,
    #[case] others: &str,
    #[case] expected: &str,
  ) {
    let tree = three_tip_tree();
    let other = Seq::try_from_str(others).unwrap();
    let sequences = btreemap! {
      tree.root  => Seq::try_from_str(root).unwrap(),
      tree.inner => other.clone(),
      tree.tip_c => other.clone(),
      tree.tip_a => other.clone(),
      tree.tip_b => other,
    };
    let result = realized(&tree.graph, &sequences);
    assert_error!(result, expected);
  }

  #[test]
  fn test_branch_lengths_perturbed_without_noise_keeps_lengths() {
    let tree = three_tip_tree();
    let lengths = btreemap! { tree.root_inner => 0.1, tree.root_c => 0.2, tree.inner_a => 0.3, tree.inner_b => 0.4 };
    let actual = perturbed(&lengths, 0.0, &mut get_random_number_generator(8)).unwrap();
    assert_eq!(lengths, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::negative(-0.1,          "Perturbation sigma must be finite and non-negative, got -0.1")]
  #[case::nan(     f64::NAN,      "Perturbation sigma must be finite and non-negative, got NaN")]
  #[case::infinite(f64::INFINITY, "Perturbation sigma must be finite and non-negative, got inf")]
  #[trace]
  fn test_branch_lengths_perturbed_rejects_invalid_sigma(#[case] sigma: f64, #[case] expected: &str) {
    let tree = three_tip_tree();
    let lengths = btreemap! { tree.root_inner => 0.1 };
    let result = perturbed(&lengths, sigma, &mut get_random_number_generator(0));
    assert_error!(result, expected);
  }
}
