#[cfg(test)]
mod tests {
  use crate::commands::homoplasy::result::{DrmAnnotation, MultiplicityRow, SiteBranchesRow, TaxonResult};
  use crate::results::homoplasy::{
    AMBIGUOUS_SITES_SHOWN, AmbiguousSite, HomoplasySite, HomoplasyStatistics, RecurrentMutation, SiteSubstitution,
    homoplasy_statistics,
  };
  use helpers::{drm, poisson_expected, site_hits, stats};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::o;
  use treetime_utils::{pretty_assert_abs_diff_eq, vec_of_owned};

  #[test]
  fn test_homoplasy_results_of_a_statistics_file() {
    let stats = stats(false);

    let actual = homoplasy_statistics(&stats);

    let expected = HomoplasyStatistics {
      drm_annotated: true,
      zero_based: false,
      genome_length: 100,
      total_branch_length: 0.5,
      substitutions: 7,
      distinct_substitutions: 4,
      recurrent_substitutions: 2,
      sites_hit_more_than_once: 2,
      expected_sites_hit_more_than_once: actual.expected_sites_hit_more_than_once,
      log_likelihood_difference: -3.5,
      samples_with_homoplasies: 1,
      recurrent_drm_substitutions: Some(1),
      ambiguous_changes: 4,
      indels: 3,
      site_hits: site_hits(),
      multiplicities: vec![
        MultiplicityRow {
          branches: 1,
          mutations: 2,
        },
        MultiplicityRow {
          branches: 2,
          mutations: 1,
        },
        MultiplicityRow {
          branches: 3,
          mutations: 1,
        },
      ],
      recurrent: vec![
        RecurrentMutation {
          mutation: o!("G10A"),
          position: 10,
          display_position: 10,
          branches: 3,
          terminal_branches: 2,
          branch_names: vec_of_owned!["A", "B", "node_1"],
          drm: Some(drm()),
        },
        RecurrentMutation {
          mutation: o!("T20C"),
          position: 20,
          display_position: 20,
          branches: 2,
          terminal_branches: 0,
          branch_names: vec_of_owned!["C", "D"],
          drm: None,
        },
      ],
      sites: vec![
        HomoplasySite {
          position: 10,
          display_position: 10,
          branches: 4,
          substitutions: vec![
            SiteSubstitution {
              mutation: o!("G10A"),
              branches: 3,
              branch_names: vec_of_owned!["A", "B", "node_1"],
              drm: Some(drm()),
            },
            SiteSubstitution {
              mutation: o!("G10T"),
              branches: 1,
              branch_names: vec_of_owned!["E"],
              drm: None,
            },
          ],
        },
        HomoplasySite {
          position: 20,
          display_position: 20,
          branches: 2,
          substitutions: vec![SiteSubstitution {
            mutation: o!("T20C"),
            branches: 2,
            branch_names: vec_of_owned!["C", "D"],
            drm: None,
          }],
        },
      ],
      recurrent_indels: vec![RecurrentMutation {
        mutation: o!("del:3-4:AC"),
        position: 3,
        display_position: 3,
        branches: 2,
        terminal_branches: 1,
        branch_names: vec_of_owned!["A", "C"],
        drm: None,
      }],
      taxa: stats.taxa,
      ambiguous_sites: vec![
        AmbiguousSite {
          position: 5,
          display_position: 5,
          branches: 3,
        },
        AmbiguousSite {
          position: 6,
          display_position: 6,
          branches: 1,
        },
      ],
      ambiguous_site_count: 2,
    };
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_based( false, ((10, 10), (20, 20), (3, 3), (5, 5)))]
  #[case::zero_based(true,  ((11, 10), (21, 20), (4, 3), (6, 5)))]
  #[trace]
  fn test_homoplasy_results_count_tree_positions_from_one(
    #[case] zero_based: bool,
    #[case] expected: ((usize, usize), (usize, usize), (usize, usize), (usize, usize)),
  ) {
    let results = homoplasy_statistics(&stats(zero_based));

    assert_eq!(
      expected,
      (
        (results.recurrent[0].position, results.recurrent[0].display_position),
        (results.sites[1].position, results.sites[1].display_position),
        (results.recurrent_indels[0].position, results.recurrent_indels[0].display_position),
        (results.ambiguous_sites[0].position, results.ambiguous_sites[0].display_position),
      )
    );
  }

  #[test]
  fn test_homoplasy_results_expect_sites_hit_more_than_once_from_the_whole_poisson_tail() {
    let results = homoplasy_statistics(&stats(false));

    let rate: f64 = 0.07;
    let expected = 100.0 * (1.0 - (-rate).exp() * (1.0 + rate));
    pretty_assert_abs_diff_eq!(expected, results.expected_sites_hit_more_than_once, epsilon = 1e-13);
    assert!(results.expected_sites_hit_more_than_once > (2..=4).map(poisson_expected).sum::<f64>());
  }

  #[test]
  fn test_homoplasy_results_list_one_site_row_per_site_hit_more_than_once() {
    let results = homoplasy_statistics(&stats(false));

    assert_eq!(
      (results.sites_hit_more_than_once, vec![4, 2]),
      (
        results.sites.len(),
        results.sites.iter().map(|site| site.branches).collect::<Vec<_>>()
      )
    );
  }

  #[test]
  fn test_homoplasy_results_show_the_first_ambiguous_sites_and_count_all() {
    let mut stats = stats(false);
    stats.ambiguous.sites = (1..=AMBIGUOUS_SITES_SHOWN + 50)
      .rev()
      .map(|position| SiteBranchesRow {
        position,
        branches: position,
      })
      .collect();

    let results = homoplasy_statistics(&stats);

    assert_eq!(
      (AMBIGUOUS_SITES_SHOWN, AMBIGUOUS_SITES_SHOWN + 50, Some(150), Some(51)),
      (
        results.ambiguous_sites.len(),
        results.ambiguous_site_count,
        results.ambiguous_sites.first().map(|site| site.branches),
        results.ambiguous_sites.last().map(|site| site.branches),
      )
    );
  }

  #[test]
  fn test_homoplasy_results_without_drms_count_no_drm_substitutions() {
    let mut stats = stats(false);
    stats.drm_annotated = false;
    for row in &mut stats.substitutions.all.ranked {
      row.drm = None;
    }

    let results = homoplasy_statistics(&stats);

    assert_eq!(
      (false, None, vec![None, None]),
      (
        results.drm_annotated,
        results.recurrent_drm_substitutions,
        results
          .recurrent
          .iter()
          .map(|row| row.drm.clone())
          .collect::<Vec<Option<DrmAnnotation>>>()
      )
    );
  }

  #[test]
  fn test_homoplasy_results_count_samples_with_a_homoplasic_substitution() {
    let mut stats = stats(false);
    stats.taxa.push(TaxonResult {
      name: o!("F"),
      homoplasic_mutations: vec_of_owned!["T20C"],
      drm_mutations: Some(0),
      ambiguous_changes: 0,
      recurrent_indels: 0,
    });

    let results = homoplasy_statistics(&stats);

    assert_eq!(2, results.samples_with_homoplasies);
  }

  mod helpers {
    use crate::commands::homoplasy::result::{
      DrmAnnotation, IndelResult, MultiplicityRow, MutationTable, RankedMutation, SiteBranchesRow, SiteHitsRow,
      SubstitutionResult, TaxonResult,
    };
    use crate::results::homoplasy::{AmbiguousCount, AmbiguousStats, HomoplasyStatsFile};
    use treetime::o;
    use treetime_utils::vec_of_owned;

    const GENOME_LENGTH: f64 = 100.0;
    const RATE: f64 = 0.07;

    pub(super) fn stats(zero_based: bool) -> HomoplasyStatsFile {
      HomoplasyStatsFile {
        zero_based,
        drm_annotated: true,
        substitutions: SubstitutionResult {
          genome_length: 100,
          total_branch_length: 0.5,
          terminal_branch_length: 0.3,
          all: table(
            7,
            &[(1, 2), (2, 1), (3, 1)],
            vec![
              ranked("G10A", 10, &["A", "B", "node_1"], Some(drm())),
              ranked("T20C", 20, &["C", "D"], None),
              ranked("G10T", 10, &["E"], None),
              ranked("A30G", 30, &["F"], None),
            ],
          ),
          terminal: table(
            3,
            &[(1, 1), (2, 1)],
            vec![
              ranked("G10A", 10, &["A", "B"], Some(drm())),
              ranked("A30G", 30, &["F"], None),
            ],
          ),
          site_hits: site_hits(),
          log_likelihood_difference: -3.5,
        },
        ambiguous: AmbiguousStats {
          all: AmbiguousCount { mutations: 4 },
          sites: vec![
            SiteBranchesRow {
              position: 5,
              branches: 3,
            },
            SiteBranchesRow {
              position: 6,
              branches: 1,
            },
          ],
        },
        indels: IndelResult {
          all: table(
            3,
            &[(1, 1), (2, 1)],
            vec![
              ranked("del:3-4:AC", 3, &["A", "C"], None),
              ranked("ins:50-50:T", 50, &["B"], None),
            ],
          ),
          terminal: table(1, &[(1, 1)], vec![ranked("del:3-4:AC", 3, &["A"], None)]),
        },
        taxa: vec![
          TaxonResult {
            name: o!("A"),
            homoplasic_mutations: vec_of_owned!["G10A"],
            drm_mutations: Some(1),
            ambiguous_changes: 0,
            recurrent_indels: 1,
          },
          TaxonResult {
            name: o!("B"),
            homoplasic_mutations: vec![],
            drm_mutations: Some(0),
            ambiguous_changes: 1,
            recurrent_indels: 0,
          },
        ],
      }
    }

    pub(super) fn site_hits() -> Vec<SiteHitsRow> {
      [(0, 97), (1, 1), (2, 1), (3, 0), (4, 1)]
        .into_iter()
        .map(|(hits, sites)| SiteHitsRow {
          hits,
          sites,
          expected: poisson_expected(hits),
        })
        .collect()
    }

    #[expect(clippy::as_conversions, reason = "hit counts below 5 convert to f64 exactly")]
    pub(super) fn poisson_expected(hits: usize) -> f64 {
      let factorial: f64 = (1..=hits).map(|k| k as f64).product();
      GENOME_LENGTH * (-RATE).exp() * RATE.powi(hits as i32) / factorial
    }

    pub(super) fn drm() -> DrmAnnotation {
      DrmAnnotation {
        gene: o!("HA"),
        drug: o!("oseltamivir"),
        substitution: Some(o!("H275Y")),
      }
    }

    fn table(mutations: usize, multiplicities: &[(usize, usize)], ranked: Vec<RankedMutation>) -> MutationTable {
      MutationTable {
        mutations,
        multiplicities: multiplicities
          .iter()
          .map(|&(branches, mutations)| MultiplicityRow { branches, mutations })
          .collect(),
        ranked,
      }
    }

    fn ranked(mutation: &str, position: usize, branches: &[&str], drm: Option<DrmAnnotation>) -> RankedMutation {
      RankedMutation {
        mutation: o!(mutation),
        position,
        multiplicity: branches.len(),
        branches: branches.iter().map(|&name| o!(name)).collect(),
        drm,
      }
    }
  }
}
