#[cfg(test)]
mod tests {
  use crate::__tests__::test_tree_output::tests::helpers::{
    Mutations, NUC_ALPHABET, ancestral_graph, ancestral_setup, auspice, dated_graph, dated_setup, mutation_on,
  };
  use crate::annotated_graph::{AnnotatedGraph, TreeSequences};
  use crate::output_plan::CommandKind;
  use eyre::Report;
  use helpers::four_leaf_auspice_color_by;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::Seq;
  use treetime_utils::o;

  #[rustfmt::skip]
  #[rstest]
  #[case::most_branches_wins(            "A=C2T,G5A B=G5T C=G5A D=C2T,A8G", Some("gt-nuc_5"))]
  #[case::tie_takes_lowest_position(     "A=G5A,C2T C=G5A D=C2T,A8G",       Some("gt-nuc_2"))]
  #[case::one_branch_per_site(           "B=A8G C=G5A",                     Some("gt-nuc_5"))]
  #[case::ambiguous_changes_do_not_count("A=G5R B=G5R C=G5N D=A8G",         Some("gt-nuc_8"))]
  #[case::only_ambiguous_changes(        "A=G5R B=G5R",                     None)]
  #[case::no_mutations(                  "",                                None)]
  #[trace]
  fn test_auspice_color_by_defaults_to_the_site_on_most_branches(
    #[case] mutations: &str,
    #[case] expected: Option<&str>,
  ) -> Result<(), Report> {
    let actual = four_leaf_auspice_color_by(mutations)?;
    assert_eq!(expected.map(str::to_owned), actual);
    Ok(())
  }

  #[rstest]
  #[case::indel(Mutations::Indel)]
  #[case::amino_acid(Mutations::AminoAcid)]
  #[trace]
  fn test_auspice_color_by_needs_a_nucleotide_substitution(#[case] mutations: Mutations) -> Result<(), Report> {
    let setup = ancestral_setup(mutations)?;
    let actual = auspice(&ancestral_graph(&setup), CommandKind::Ancestral)?;
    assert_eq!(None, actual.data.meta.display_defaults.color_by);
    Ok(())
  }

  #[test]
  fn test_auspice_color_by_keeps_the_filter_default_over_genotype() -> Result<(), Report> {
    let setup = dated_setup()?;
    let root_sequence = Seq::try_from_str("ACGT")?;
    let edge_mutations = mutation_on(&setup.topology, "A")?;
    let graph = AnnotatedGraph {
      sequences: Some(TreeSequences {
        alphabet: &NUC_ALPHABET,
        root_sequence: &root_sequence,
        edge_mutations: &edge_mutations,
        mutation_counts: None,
        amino_acids: None,
      }),
      ..dated_graph(&setup, None)
    };
    let actual = auspice(&graph, CommandKind::Timetree)?;
    assert_eq!(Some(o!("bad_branch")), actual.data.meta.display_defaults.color_by);
    Ok(())
  }

  mod helpers {
    use crate::__tests__::test_tree_output::tests::helpers::{
      NUC_ALPHABET, annotated, auspice, no_mutations, parent_edge, topology_from,
    };
    use crate::annotated_graph::{AnnotatedGraph, TreeSequences};
    use crate::output_plan::CommandKind;
    use eyre::Report;
    use std::str::FromStr;
    use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
    use treetime_primitives::Seq;

    pub(super) fn four_leaf_auspice_color_by(mutations: &str) -> Result<Option<String>, Report> {
      let topology = topology_from("((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root;")?;
      let mut edge_mutations = no_mutations(&topology);
      for branch in mutations.split_whitespace() {
        let (name, subs) = branch.split_once('=').expect("fixture branch must be NAME=SUBS");
        let subs = subs
          .split(',')
          .map(|sub| Ok(Mutation::substitution(MutationTrack::Nucleotide, Sub::from_str(sub)?)))
          .collect::<Result<_, Report>>()?;
        edge_mutations.insert(parent_edge(&topology, name)?, subs);
      }
      let root_sequence = Seq::try_from_str("ACGTGCAAAC")?;
      let graph = AnnotatedGraph {
        sequences: Some(TreeSequences {
          alphabet: &NUC_ALPHABET,
          root_sequence: &root_sequence,
          edge_mutations: &edge_mutations,
          mutation_counts: None,
          amino_acids: None,
        }),
        ..annotated(&topology)
      };
      Ok(
        auspice(&graph, CommandKind::Ancestral)?
          .data
          .meta
          .display_defaults
          .color_by,
      )
    }
  }
}
