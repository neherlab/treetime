#[cfg(test)]
mod tests {
  use crate::commands::homoplasy::drms::{DrmRow, DrmTable};
  use crate::commands::homoplasy::result::{DrmAnnotation, RankedMutation, TaxonResult};
  use crate::commands::homoplasy::summary::{ResultContext, homoplasy_result};
  use eyre::Report;
  use helpers::{names, output};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::o;
  use treetime_utils::vec_of_owned;

  #[rustfmt::skip]
  #[rstest]
  #[case::one_based( false, ("G10A", "del:3-4:AC", 5))]
  #[case::zero_based(true,  ("G9A",  "del:2-3:AC", 4))]
  #[trace]
  fn test_summary_homoplasy_positions_follow_indexing(
    #[case] zero_based: bool,
    #[case] (substitution, indel, ambiguous_site): (&str, &str, usize),
  ) -> Result<(), Report> {
    let names = names();
    let result = homoplasy_result(&output()?, &ResultContext { names: &names, drms: None, zero_based });

    assert_eq!(
      (o!(substitution), o!(indel), ambiguous_site),
      (
        result.substitutions.all.ranked[0].mutation.clone(),
        result.indels.all.ranked[0].mutation.clone(),
        result.ambiguous.sites[0].position
      )
    );
    Ok(())
  }

  #[test]
  fn test_summary_homoplasy_names_branches_and_annotates_drms() -> Result<(), Report> {
    let names = names();
    let drms = DrmTable::from_rows(vec![DrmRow {
      genomic_position: 10,
      alt_base: o!("A"),
      drug: o!("NRTI"),
      gene: o!("RT"),
      substitution: o!("M41L"),
    }])?;
    let result = homoplasy_result(
      &output()?,
      &ResultContext {
        names: &names,
        drms: Some(&drms),
        zero_based: false,
      },
    );

    let expected = RankedMutation {
      mutation: o!("G10A"),
      multiplicity: 2,
      branches: vec_of_owned!["A", "node_2"],
      drm: Some(DrmAnnotation {
        gene: o!("RT"),
        drug: o!("NRTI"),
        substitution: Some(o!("M41L")),
      }),
    };
    assert_eq!(expected, result.substitutions.all.ranked[0]);
    assert_eq!(None, result.ambiguous.all.ranked[0].drm);
    Ok(())
  }

  #[test]
  fn test_summary_homoplasy_orders_taxa_by_homoplasic_count_then_name() -> Result<(), Report> {
    let names = names();
    let drms = DrmTable::from_rows(vec![DrmRow {
      genomic_position: 10,
      alt_base: o!("T"),
      drug: o!("NRTI"),
      gene: o!("RT"),
      substitution: o!("M41L"),
    }])?;
    let result = homoplasy_result(
      &output()?,
      &ResultContext {
        names: &names,
        drms: Some(&drms),
        zero_based: false,
      },
    );

    let expected = vec![
      TaxonResult {
        name: o!("A"),
        homoplasic_mutations: vec_of_owned!["G10A"],
        drm_mutations: Some(1),
        ambiguous_changes: 0,
        recurrent_indels: 1,
      },
      TaxonResult {
        name: o!("B"),
        homoplasic_mutations: vec_of_owned!["C20T"],
        drm_mutations: Some(0),
        ambiguous_changes: 2,
        recurrent_indels: 0,
      },
      TaxonResult {
        name: o!("C"),
        homoplasic_mutations: vec![],
        drm_mutations: Some(0),
        ambiguous_changes: 1,
        recurrent_indels: 0,
      },
    ];
    assert_eq!(expected, result.taxa);
    Ok(())
  }

  mod helpers {
    use eyre::Report;
    use maplit::btreemap;
    use std::collections::BTreeMap;
    use std::str::FromStr;
    use treetime::homoplasy::pipeline::{
      AmbiguousStats, HomoplasyOutput, IndelKey, IndelStats, SiteBranches, SubstitutionStats,
    };
    use treetime::homoplasy::recurrence::{Recurrence, RecurrenceTable};
    use treetime::homoplasy::site_hits::{SiteHistogram, SiteHits};
    use treetime::o;
    use treetime::seq::indel::InDelKind;
    use treetime::seq::mutation::Sub;
    use treetime_graph::node::GraphNodeKey;
    use treetime_primitives::Seq;

    const A: GraphNodeKey = GraphNodeKey(1);
    const UNNAMED: GraphNodeKey = GraphNodeKey(2);
    const B: GraphNodeKey = GraphNodeKey(3);
    const C: GraphNodeKey = GraphNodeKey(4);

    pub(super) fn names() -> BTreeMap<GraphNodeKey, Option<String>> {
      btreemap! {
        A => Some(o!("A")),
        UNNAMED => None,
        B => Some(o!("B")),
        C => Some(o!("C")),
      }
    }

    pub(super) fn output() -> Result<HomoplasyOutput, Report> {
      let recurrent = IndelKey {
        kind: InDelKind::Deletion,
        range: (2, 4),
        sequence: Seq::try_from_slice(b"AC")?,
      };
      Ok(HomoplasyOutput {
        substitutions: SubstitutionStats {
          total_branch_length: 1.0,
          terminal_branch_length: 0.5,
          all: table(vec![(Sub::from_str("G10A")?, vec![A, UNNAMED]), (Sub::from_str("C20T")?, vec![B])]),
          terminal: table(vec![(Sub::from_str("G10A")?, vec![A]), (Sub::from_str("C20T")?, vec![B])]),
          sites: SiteHistogram {
            genome_length: 30,
            rows: vec![SiteHits {
              hits: 0,
              sites: 28,
              expected: 28.0,
            }],
            log_likelihood_difference: 0.0,
          },
          leaves: btreemap! {
            A => vec![Sub::from_str("G10A")?],
            B => vec![Sub::from_str("C20T")?],
          },
        },
        ambiguous: AmbiguousStats {
          all: table(vec![(Sub::from_str("G5R")?, vec![B, C])]),
          sites: vec![SiteBranches {
            position: 4,
            branches: 2,
          }],
          leaves: btreemap! {B => 2, C => 1},
        },
        indels: IndelStats {
          all: table(vec![(recurrent.clone(), vec![A, UNNAMED])]),
          terminal: table(vec![(recurrent, vec![A])]),
          leaves: btreemap! {A => 1},
        },
      })
    }

    fn table<K>(ranked: Vec<(K, Vec<GraphNodeKey>)>) -> RecurrenceTable<K> {
      let ranked: Vec<Recurrence<K>> = ranked
        .into_iter()
        .map(|(mutation, branches)| Recurrence { mutation, branches })
        .collect();
      let mut histogram = BTreeMap::new();
      for recurrence in &ranked {
        *histogram.entry(recurrence.multiplicity()).or_default() += 1;
      }
      RecurrenceTable {
        count: ranked.iter().map(Recurrence::multiplicity).sum(),
        histogram,
        ranked,
      }
    }
  }
}
