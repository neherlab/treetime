#[cfg(test)]
mod tests {
  use crate::config::catalog::SettingPurpose;
  use crate::json_float::JsonFloat;
  use crate::results::clades::{clade_in_trees, matched_ancestors, same_topology};
  use crate::results::compare::{
    AncestorComparison, AncestorShift, LikelihoodMismatch, compare_ancestors, compare_estimates, likelihood_mismatches,
  };
  use crate::results::timetree::{CoalescentPrior, TimetreeEstimates};
  use crate::results::tree::DateInterval;
  use crate::results::year_date::YearDate;
  use crate::runs::setting_differences::{DifferenceValues, SettingDifference};
  use helpers::{difference, estimates, resolved, tree};
  use pretty_assertions::assert_eq;
  use treetime::o;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_clades_match_across_tip_order() {
    let first = tree(&[
      ("root", None, Some(2000.0)),
      ("X", Some(0), Some(2005.0)),
      ("A", Some(1), None),
      ("B", Some(1), None),
      ("C", Some(0), None),
    ]);
    let second = tree(&[
      ("root", None, Some(2001.0)),
      ("C", Some(0), None),
      ("Y", Some(0), Some(2005.5)),
      ("B", Some(2), None),
      ("A", Some(2), None),
    ]);

    assert_eq!(vec![(0, 0), (1, 2)], matched_ancestors(&first, &second));
  }

  #[test]
  fn test_clades_present_in_one_tree_only_do_not_match() {
    let first = tree(&[
      ("root", None, None),
      ("X", Some(0), None),
      ("A", Some(1), None),
      ("B", Some(1), None),
      ("C", Some(0), None),
    ]);
    let other = tree(&[
      ("root", None, None),
      ("Z", Some(0), None),
      ("A", Some(1), None),
      ("C", Some(1), None),
      ("B", Some(0), None),
    ]);

    assert_eq!(vec![(0, 0)], matched_ancestors(&first, &other));
    assert_eq!(vec![None], clade_in_trees(&first, 1, [&other]));
  }

  #[test]
  fn test_clade_in_trees_finds_a_clade_and_a_sample_in_each_tree() {
    let first = tree(&[
      ("root", None, None),
      ("X", Some(0), None),
      ("A", Some(1), None),
      ("B", Some(1), None),
      ("C", Some(0), None),
    ]);
    let reordered = tree(&[
      ("root", None, None),
      ("C", Some(0), None),
      ("Y", Some(0), None),
      ("B", Some(2), None),
      ("A", Some(2), None),
    ]);
    let other = tree(&[
      ("root", None, None),
      ("Z", Some(0), None),
      ("A", Some(1), None),
      ("C", Some(1), None),
      ("B", Some(0), None),
    ]);

    assert_eq!(
      (vec![Some(2), None], vec![Some(4), Some(2)]),
      (
        clade_in_trees(&first, 1, [&reordered, &other]),
        clade_in_trees(&first, 2, [&reordered, &other])
      )
    );
  }

  #[test]
  fn test_compare_ancestors_shifts_dates_of_shared_ancestors_in_days() {
    let first = tree(&[
      ("root", None, Some(2000.0)),
      ("X", Some(0), Some(2005.0)),
      ("A", Some(1), None),
      ("B", Some(1), None),
      ("C", Some(0), None),
    ]);
    let second = tree(&[
      ("root", None, Some(2001.0)),
      ("C", Some(0), None),
      ("Y", Some(0), Some(2005.5)),
      ("B", Some(2), None),
      ("A", Some(2), None),
    ]);

    let expected = AncestorComparison {
      shifts: vec![
        AncestorShift {
          name: o!("root"),
          tips: 3,
          date_first: YearDate::new(2000.0),
          shift_days: 366.0,
        },
        AncestorShift {
          name: o!("X"),
          tips: 2,
          date_first: YearDate::new(2005.0),
          shift_days: 182.5,
        },
      ],
      ancestors: 2,
      mean_absolute_shift_days: Some(274.25),
    };
    assert_eq!(expected, compare_ancestors(&first, &second));
  }

  #[test]
  fn test_compare_estimates_reports_differences_second_minus_first() {
    let first = estimates(2000.0, (1999.0, 2001.0), 1e-3, 1, JsonFloat(-100.0));
    let second = estimates(2001.0, (2000.0, 2001.0), 1.5e-3, 3, JsonFloat(-90.0));

    let comparison = compare_estimates(first, second, vec![]);

    assert_eq!(Some(366.0), comparison.root_shift_days);
    assert_eq!(Some(366.0 - 731.0), comparison.root_interval_change_days);
    pretty_assert_ulps_eq!(50.0, comparison.clock_rate_change_percent.unwrap(), max_ulps = 4);
    assert_eq!(2, comparison.excluded_samples_change);
    assert_eq!(Some(10.0), comparison.log_likelihood_change);
  }

  #[test]
  fn test_compare_estimates_leaves_out_a_non_finite_log_likelihood() {
    let first = estimates(2000.0, (1999.0, 2001.0), 1e-3, 1, JsonFloat(f64::INFINITY));
    let second = estimates(2001.0, (2000.0, 2001.0), 1e-3, 1, JsonFloat(-90.0));

    assert_eq!(None, compare_estimates(first, second, vec![]).log_likelihood_change);
  }

  #[test]
  fn test_likelihood_mismatches_of_the_same_model_and_tree_are_empty() {
    let first = estimates(2000.0, (1999.0, 2001.0), 1e-3, 1, JsonFloat(-100.0));
    let second = estimates(2001.0, (2000.0, 2001.0), 2e-3, 1, JsonFloat(-90.0));
    let mismatches = likelihood_mismatches(&first, &second, &[], Some((&resolved(), &resolved())));

    assert_eq!(
      (vec![], Some(10.0)),
      (
        mismatches.clone(),
        compare_estimates(first, second, mismatches).log_likelihood_change
      )
    );
  }

  #[test]
  fn test_likelihood_mismatches_of_another_polytomy_resolution_name_the_topology() {
    let first = estimates(2000.0, (1999.0, 2001.0), 1e-3, 1, JsonFloat(-100.0));
    let second = estimates(2000.0, (1999.0, 2001.0), 1e-3, 1, JsonFloat(-102.0));
    let other = tree(&[
      ("root", None, None),
      ("Y", Some(0), None),
      ("A", Some(1), None),
      ("C", Some(1), None),
      ("B", Some(0), None),
    ]);
    let mismatches = likelihood_mismatches(&first, &second, &[], Some((&resolved(), &other)));

    assert_eq!(
      (vec![LikelihoodMismatch::Topology], None),
      (
        mismatches.clone(),
        compare_estimates(first, second, mismatches).log_likelihood_change
      )
    );
  }

  #[test]
  fn test_likelihood_mismatches_of_runs_without_trees_name_the_topology() {
    let first = estimates(2000.0, (1999.0, 2001.0), 1e-3, 1, JsonFloat(-100.0));
    assert_eq!(
      vec![LikelihoodMismatch::Topology],
      likelihood_mismatches(&first, &first.clone(), &[], None)
    );
  }

  #[test]
  fn test_likelihood_mismatches_name_the_coalescent_prior_the_model_and_the_data() {
    let first = estimates(2000.0, (1999.0, 2001.0), 1e-3, 1, JsonFloat(-100.0));
    let second = TimetreeEstimates {
      coalescent_prior: CoalescentPrior::Optimized,
      ..first.clone()
    };
    let differences = [
      difference("tree", "Input data", SettingPurpose::Data),
      difference("model", "Substitution model", SettingPurpose::Analysis),
      difference("clock_rate", "Molecular clock", SettingPurpose::Analysis),
    ];

    assert_eq!(
      vec![
        LikelihoodMismatch::Data,
        LikelihoodMismatch::SubstitutionModel,
        LikelihoodMismatch::CoalescentPrior,
      ],
      likelihood_mismatches(&first, &second, &differences, Some((&resolved(), &resolved())))
    );
  }

  #[test]
  fn test_same_topology_ignores_the_order_of_children() {
    let reordered = tree(&[
      ("root", None, None),
      ("C", Some(0), None),
      ("X", Some(0), None),
      ("B", Some(2), None),
      ("A", Some(2), None),
    ]);
    assert_eq!(
      (true, false),
      (
        same_topology(&resolved(), &reordered),
        same_topology(&resolved(), &tree(&[("root", None, None), ("A", Some(0), None)]))
      )
    );
  }

  mod helpers {
    use super::*;
    use crate::results::tree::{ResultNode, ResultTree};
    use crate::results::year_date::YearDate;

    pub(super) fn tree(spec: &[(&str, Option<usize>, Option<f64>)]) -> ResultTree {
      let mut nodes = spec
        .iter()
        .map(|&(name, parent, date)| ResultNode {
          name: name.to_owned(),
          parent,
          children: vec![],
          tips: 0,
          div: None,
          date: date.map(YearDate::new),
          date_interval: None,
          excluded: None,
          mutations: vec![],
        })
        .collect::<Vec<_>>();
      for index in 1..nodes.len() {
        let parent = nodes[index].parent.unwrap();
        nodes[parent].children.push(index);
      }
      for index in (0..nodes.len()).rev() {
        nodes[index].tips = if nodes[index].children.is_empty() {
          1
        } else {
          nodes[index].children.iter().map(|&child| nodes[child].tips).sum()
        };
      }
      ResultTree {
        nodes,
        colorings: vec![],
        default_color_by: None,
      }
    }

    pub(super) fn resolved() -> ResultTree {
      tree(&[
        ("root", None, None),
        ("X", Some(0), None),
        ("A", Some(1), None),
        ("B", Some(1), None),
        ("C", Some(0), None),
      ])
    }

    pub(super) fn difference(key: &str, group: &str, purpose: SettingPurpose) -> SettingDifference {
      SettingDifference {
        key: key.to_owned(),
        label: key.to_owned(),
        group: group.to_owned(),
        purpose,
        values: DifferenceValues::Setting {
          first: None,
          second: None,
        },
      }
    }

    pub(super) fn estimates(
      root_date: f64,
      (lower, upper): (f64, f64),
      clock_rate: f64,
      excluded_samples: usize,
      log_likelihood: JsonFloat,
    ) -> TimetreeEstimates {
      TimetreeEstimates {
        root_date: Some(YearDate::new(root_date)),
        root_interval: DateInterval::new(lower, upper),
        root_near_interval_edge: false,
        clock_rate: Some(clock_rate),
        clock_rate_std: None,
        clock_rate_fixed: false,
        r: None,
        r_squared: None,
        samples: 10,
        latest_sample: None,
        excluded_samples,
        coalescent_prior: CoalescentPrior::None,
        relaxed_clock: None,
        log_likelihood: Some(log_likelihood),
        iterations: 2,
      }
    }
  }
}
