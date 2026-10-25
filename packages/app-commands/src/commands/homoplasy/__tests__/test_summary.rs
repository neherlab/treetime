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
  #[case::one_based( false, (("G10A", 10), ("del:3-4:AC", 3), 5))]
  #[case::zero_based(true,  (("G9A",  9),  ("del:2-3:AC", 2), 4))]
  #[trace]
  fn test_summary_homoplasy_positions_follow_indexing(
    #[case] zero_based: bool,
    #[case] ((substitution, substitution_position), (indel, indel_position), ambiguous_site): ((&str, usize), (&str, usize), usize),
  ) -> Result<(), Report> {
    let (names, taxa) = names()?;
    let result = homoplasy_result(&output(taxa)?, &ResultContext { names: &names, drms: None, zero_based });
    let first_substitution = &result.substitutions.all.ranked[0];
    let first_indel = &result.indels.all.ranked[0];

    assert_eq!(
      ((o!(substitution), substitution_position), (o!(indel), indel_position), ambiguous_site),
      (
        (first_substitution.mutation.clone(), first_substitution.position),
        (first_indel.mutation.clone(), first_indel.position),
        result.ambiguous.sites[0].position
      )
    );
    Ok(())
  }

  #[test]
  fn test_summary_homoplasy_names_branches_and_annotates_drms() -> Result<(), Report> {
    let (names, taxa) = names()?;
    let drms = DrmTable::from_rows(vec![DrmRow {
      genomic_position: 10,
      alt_base: o!("A"),
      drug: o!("NRTI"),
      gene: o!("RT"),
      substitution: o!("M41L"),
    }])?;
    let result = homoplasy_result(
      &output(taxa)?,
      &ResultContext {
        names: &names,
        drms: Some(&drms),
        zero_based: false,
      },
    );

    let expected = RankedMutation {
      mutation: o!("G10A"),
      position: 10,
      multiplicity: 2,
      branches: vec_of_owned!["A", "NODE_0000000"],
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
    let (names, taxa) = names()?;
    let drms = DrmTable::from_rows(vec![DrmRow {
      genomic_position: 10,
      alt_base: o!("T"),
      drug: o!("NRTI"),
      gene: o!("RT"),
      substitution: o!("M41L"),
    }])?;
    let result = homoplasy_result(
      &output(taxa)?,
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
    use treetime_graph::keymap;

    use std::collections::BTreeMap;
    use std::str::FromStr;
    use treetime::homoplasy::pipeline::{AmbiguousStats, HomoplasyOutput, IndelKey, IndelStats, SubstitutionStats};
    use treetime::homoplasy::recurrence::{Recurrence, RecurrenceTable};
    use treetime::homoplasy::site_branches::SiteBranches;
    use treetime::homoplasy::site_hits::{SiteHistogram, SiteHits};
    use treetime::o;
    use treetime::seq::indel::InDelKind;
    use treetime::seq::mutation::Sub;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_primitives::Seq;

    #[derive(Clone, Copy)]
    pub(super) struct Taxa {
      a: GraphNodeKey,
      unnamed: GraphNodeKey,
      b: GraphNodeKey,
      c: GraphNodeKey,
    }

    pub(super) fn names() -> Result<(NodeNames, Taxa), Report> {
      let mut graph = Graph::new();
      let root = graph.add_node();
      let [a, unnamed, b, c] = [(); 4].map(|()| graph.add_node());
      for child in [a, unnamed, b, c] {
        graph.add_edge(root, child)?;
      }
      graph.build()?;
      let input = keymap! {
        root => Some(o!("root")),
        a => Some(o!("A")),
        unnamed => None,
        b => Some(o!("B")),
        c => Some(o!("C")),
      };
      Ok((NodeNames::from_input(input, &graph)?.names, Taxa { a, unnamed, b, c }))
    }

    pub(super) fn output(taxa: Taxa) -> Result<HomoplasyOutput, Report> {
      let Taxa { a, unnamed, b, c } = taxa;
      let recurrent = IndelKey {
        kind: InDelKind::Deletion,
        range: (2, 4),
        sequence: Seq::try_from_slice(b"AC")?,
      };
      Ok(HomoplasyOutput {
        substitutions: SubstitutionStats {
          total_branch_length: 1.0,
          terminal_branch_length: 0.5,
          all: table(vec![
            (Sub::from_str("G10A")?, vec![a, unnamed]),
            (Sub::from_str("C20T")?, vec![b]),
          ]),
          terminal: table(vec![
            (Sub::from_str("G10A")?, vec![a]),
            (Sub::from_str("C20T")?, vec![b]),
          ]),
          sites: SiteHistogram {
            genome_length: 30,
            rows: vec![SiteHits {
              hits: 0,
              sites: 28,
              expected: 28.0,
            }],
            log_likelihood_difference: 0.0,
          },
          leaves: keymap! {
            a => vec![Sub::from_str("G10A")?],
            b => vec![Sub::from_str("C20T")?],
          },
        },
        ambiguous: AmbiguousStats {
          all: table(vec![(Sub::from_str("G5R")?, vec![b, c])]),
          sites: vec![SiteBranches {
            position: 4,
            branches: 2,
          }],
          leaves: keymap! {b => 2, c => 1},
        },
        indels: IndelStats {
          all: table(vec![(recurrent.clone(), vec![a, unnamed])]),
          terminal: table(vec![(recurrent, vec![a])]),
          leaves: keymap! {a => 1},
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
