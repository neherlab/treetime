#[cfg(test)]
mod tests {
  use crate::commands::homoplasy::report::render_homoplasy_report;
  use crate::commands::homoplasy::result::DrmAnnotation;
  use helpers::fixture;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime::o;

  #[test]
  fn test_report_homoplasy_default_layout() {
    let expected = indoc! {r#"
    The TOTAL tree length is 6.617e-02 and 7 mutations were observed.
    Of these 7 mutations,
    	 - 2 occur 1 times
    	 - 1 occur 2 times
    	 - 1 occur 3 times

    Of the 100 positions in the genome,
    	 - 96 were hit 0 times (expected 93.24)
    	 - 2 were hit 1 times (expected 6.53)
    	 - 1 were hit 2 times (expected 0.23)
    	 - 1 were hit 3 times (expected 0.01)

    log-likelihood difference to Poisson distribution with same mean: -7.077e+00


    The 10 most homoplasic mutations are:
    	mut	multiplicity
    	G10A	3
    	T20C	2


    Changes involving ambiguous characters: 3 were observed.
    Of these 3 changes,
    	 - 1 occur 1 times
    	 - 1 occur 2 times

    The 10 most frequent changes involving ambiguous characters are:
    	mut	multiplicity
    	G5R	2


    Insertions and deletions: 2 were observed.
    Of these 2 insertions and deletions,
    	 - 1 occur 2 times

    The 10 most frequent insertions and deletions are:
    	indel	multiplicity
    	del:3-4:AC	2
"#};
    assert_eq!(expected, render_homoplasy_report(&fixture(), 10, false));
  }

  #[test]
  fn test_report_homoplasy_detailed_with_drms_and_one_row() {
    let mut result = fixture();
    result.drm_annotated = true;
    let drm = DrmAnnotation {
      gene: o!("RT"),
      drug: o!("NRTI"),
      substitution: Some(o!("M41L")),
    };
    result.substitutions.all.ranked[0].drm = Some(drm.clone());
    result.substitutions.terminal.ranked[0].drm = Some(drm);
    result.taxa[0].drm_mutations = Some(1);

    let expected = indoc! {r#"
    The TOTAL tree length is 6.617e-02 and 7 mutations were observed.
    Of these 7 mutations,
    	 - 2 occur 1 times
    	 - 1 occur 2 times
    	 - 1 occur 3 times

    The TERMINAL branch length is 3.962e-02 and 3 mutations were observed.
    Of these 3 mutations,
    	 - 1 occur 1 times
    	 - 1 occur 2 times

    Of the 100 positions in the genome,
    	 - 96 were hit 0 times (expected 93.24)
    	 - 2 were hit 1 times (expected 6.53)
    	 - 1 were hit 2 times (expected 0.23)
    	 - 1 were hit 3 times (expected 0.01)

    log-likelihood difference to Poisson distribution with same mean: -7.077e+00


    The 1 most homoplasic mutations are:
    	mut	multiplicity	DRM details (gene drug AAmut)
    	G10A	3	RT NRTI M41L


    The 1 most homoplasic mutations on terminal branches are:
    	mut	multiplicity	DRM details (gene drug AAmut)
    	G10A	2	RT NRTI M41L


    Taxons that carry positions that mutated elsewhere in the tree:
    	taxon name	#of homoplasic mutations	# DRM	# ambiguous changes	# recurrent indels
    	C	2	1	1	1


    Changes involving ambiguous characters: 3 were observed.
    Of these 3 changes,
    	 - 1 occur 1 times
    	 - 1 occur 2 times

    The 1 most frequent changes involving ambiguous characters are:
    	mut	multiplicity
    	G5R	2


    Insertions and deletions: 2 were observed.
    Of these 2 insertions and deletions,
    	 - 1 occur 2 times

    The 1 most frequent insertions and deletions are:
    	indel	multiplicity
    	del:3-4:AC	2
"#};
    assert_eq!(expected, render_homoplasy_report(&result, 1, true));
  }

  #[test]
  fn test_report_homoplasy_drm_without_listed_substitution() {
    let mut result = fixture();
    result.drm_annotated = true;
    result.substitutions.all.ranked[0].drm = Some(DrmAnnotation {
      gene: o!("RT"),
      drug: o!("NRTI"),
      substitution: None,
    });

    let expected = indoc! {r#"
    The 10 most homoplasic mutations are:
    	mut	multiplicity	DRM details (gene drug AAmut)
    	G10A	3	RT NRTI
    	T20C	2	
"#};
    let report = render_homoplasy_report(&result, 10, false);
    let block = report
      .split("\n\n\n")
      .find(|block| block.starts_with("The 10 most homoplasic"))
      .unwrap();
    assert_eq!(expected, format!("{block}\n"));
  }

  mod helpers {
    use crate::commands::homoplasy::result::{
      AmbiguousResult, HomoplasyResult, IndelResult, MultiplicityRow, MutationTable, RankedMutation, SiteBranchesRow,
      SiteHitsRow, SubstitutionResult, TaxonResult,
    };
    use treetime::o;
    use treetime_utils::vec_of_owned;

    pub(super) fn fixture() -> HomoplasyResult {
      HomoplasyResult {
        zero_based: false,
        drm_annotated: false,
        substitutions: SubstitutionResult {
          genome_length: 100,
          total_branch_length: 0.066_17,
          terminal_branch_length: 0.039_62,
          all: table(
            7,
            &[(1, 2), (2, 1), (3, 1)],
            &[("G10A", 3), ("T20C", 2), ("A30G", 1), ("C40T", 1)],
          ),
          terminal: table(3, &[(1, 1), (2, 1)], &[("G10A", 2), ("C40T", 1)]),
          site_hits: vec![
            site_hits(0, 96, 93.244),
            site_hits(1, 2, 6.526),
            site_hits(2, 1, 0.228),
            site_hits(3, 1, 0.005_3),
          ],
          log_likelihood_difference: -7.077,
        },
        ambiguous: AmbiguousResult {
          all: table(3, &[(1, 1), (2, 1)], &[("G5R", 2), ("A6N", 1)]),
          sites: vec![
            SiteBranchesRow {
              position: 5,
              branches: 2,
            },
            SiteBranchesRow {
              position: 6,
              branches: 1,
            },
          ],
        },
        indels: IndelResult {
          all: table(2, &[(2, 1)], &[("del:3-4:AC", 2)]),
          terminal: table(1, &[(1, 1)], &[("del:3-4:AC", 1)]),
        },
        taxa: vec![TaxonResult {
          name: o!("C"),
          homoplasic_mutations: vec_of_owned!["G10A", "T20C"],
          drm_mutations: None,
          ambiguous_changes: 1,
          recurrent_indels: 1,
        }],
      }
    }

    fn table(mutations: usize, multiplicities: &[(usize, usize)], ranked: &[(&str, usize)]) -> MutationTable {
      MutationTable {
        mutations,
        multiplicities: multiplicities
          .iter()
          .map(|&(branches, mutations)| MultiplicityRow { branches, mutations })
          .collect(),
        ranked: ranked
          .iter()
          .map(|&(mutation, multiplicity)| RankedMutation {
            mutation: mutation.to_owned(),
            multiplicity,
            branches: vec![],
            drm: None,
          })
          .collect(),
      }
    }

    fn site_hits(hits: usize, sites: usize, expected: f64) -> SiteHitsRow {
      SiteHitsRow { hits, sites, expected }
    }
  }
}
