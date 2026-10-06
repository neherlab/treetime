#[cfg(test)]
mod tests {
  use crate::__tests__::test_tree_output::tests::helpers::{
    NUC_ALPHABET, annotated, dated_graph, dated_setup, mugration_graph, mugration_setup, no_mutations, node_key,
    parent_edge, substitution, topology_from,
  };
  use crate::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, TreeSequences};
  use crate::nwk_comments::nwk_node_comments;
  use crate::output_plan::{CommandKind, TreeWriteKind};
  use crate::tree_output::write_tree_outputs;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use tempfile::TempDir;
  use treetime::progress::NoopProgress;
  use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack};
  use treetime_io::nwk::{NewickValue, NwkStyle};
  use treetime_primitives::Seq;
  use treetime_utils::io::fs::read_file_to_string;
  use treetime_utils::o;

  #[test]
  fn test_nwk_comments_mutations_are_1_based_with_indels_expanded() -> Result<(), Report> {
    let mutations = vec![
      substitution(b'A', 0, b'T')?,
      substitution(b'G', 5, b'C')?,
      Mutation {
        track: MutationTrack::Nucleotide,
        event: MutationEvent::Deletion(AlignedMutation {
          range: (1, 3),
          sequence: Seq::try_from_str("CG")?,
        }),
      },
    ];

    let actual = helpers::leaf_comments(mutations)?;

    assert_eq!(
      vec![(o!("mutations"), NewickValue::String(o!("A1T,C2-,G3-,G6C")))],
      actual
    );
    Ok(())
  }

  #[test]
  fn test_nwk_comments_mutations_sorted_by_position() -> Result<(), Report> {
    let mutations = vec![
      substitution(b'C', 50, b'G')?,
      substitution(b'A', 10, b'T')?,
      substitution(b'G', 30, b'C')?,
    ];

    let actual = helpers::leaf_comments(mutations)?;

    assert_eq!(
      vec![(o!("mutations"), NewickValue::String(o!("A11T,G31C,C51G")))],
      actual
    );
    Ok(())
  }

  #[test]
  fn test_nwk_comments_edge_without_mutations_has_none() -> Result<(), Report> {
    assert_eq!(Vec::<(String, NewickValue)>::new(), helpers::leaf_comments(vec![])?);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_root_has_no_mutation_comment() -> Result<(), Report> {
    let topology = topology_from("(A:0.1)root;")?;
    let mut edge_mutations = no_mutations(&topology);
    edge_mutations.insert(parent_edge(&topology, "A")?, vec![substitution(b'A', 0, b'T')?]);
    let root_sequence = Seq::try_from_str("ACGT")?;
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

    let comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;

    assert_eq!(
      Vec::<(String, NewickValue)>::new(),
      comments[&node_key(&topology, "root")]
    );
    Ok(())
  }

  #[test]
  fn test_nwk_comments_dates_without_sequences_give_date_comments() -> Result<(), Report> {
    let setup = dated_setup()?;
    let graph = dated_graph(&setup, None);

    let comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;

    let actual: BTreeMap<String, Vec<(String, NewickValue)>> = comments
      .into_iter()
      .map(|(key, comments)| (setup.topology.names[&key].clone().unwrap(), comments))
      .collect();
    let expected = btreemap! {
      o!("A") => vec![(o!("date"), NewickValue::NumberText(o!("2020.00")))],
      o!("B") => vec![(o!("date"), NewickValue::NumberText(o!("2021.00")))],
      o!("C") => vec![(o!("date"), NewickValue::NumberText(o!("2022.00")))],
      o!("root") => vec![(o!("date"), NewickValue::NumberText(o!("2023.00")))],
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_sequences_and_dates_give_both_comments() -> Result<(), Report> {
    let setup = dated_setup()?;
    let mut edge_mutations = no_mutations(&setup.topology);
    edge_mutations.insert(parent_edge(&setup.topology, "A")?, vec![substitution(b'A', 0, b'T')?]);
    let root_sequence = Seq::try_from_str("ACGT")?;
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

    let comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;

    let expected = vec![
      (o!("mutations"), NewickValue::String(o!("A1T"))),
      (o!("date"), NewickValue::NumberText(o!("2020.00"))),
    ];
    assert_eq!(expected, comments[&node_key(&setup.topology, "A")]);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_traits_give_trait_comments() -> Result<(), Report> {
    let setup = mugration_setup()?;
    let graph = mugration_graph(&setup, "country");

    let comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;

    let actual: BTreeMap<String, Vec<(String, NewickValue)>> = comments
      .into_iter()
      .map(|(key, comments)| (setup.topology.names[&key].clone().unwrap(), comments))
      .collect();
    let expected = btreemap! {
      o!("A") => vec![(o!("country"), NewickValue::String(o!("CH")))],
      o!("B") => vec![(o!("country"), NewickValue::String(o!("US")))],
      o!("C") => vec![(o!("country"), NewickValue::String(o!("CH")))],
      o!("root") => vec![(o!("country"), NewickValue::String(o!("US")))],
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::leading_zero("01",        r#""01""#)]
  #[case::nan_word(    "Nan",       r#""Nan""#)]
  #[case::inf_word(    "Inf",       r#""Inf""#)]
  #[case::exponent(    "1e3",       r#""1e3""#)]
  #[case::boolean_word("true",      r#""true""#)]
  #[case::comma(       "Congo, DR", r#""Congo, DR""#)]
  #[trace]
  fn test_nwk_comments_trait_values_keep_their_text(#[case] value: &str, #[case] expected: &str) -> Result<(), Report> {
    let mut setup = mugration_setup()?;
    let a_key = node_key(&setup.topology, "A");
    setup.values.insert(a_key, Some(o!(value)));
    let graph = mugration_graph(&setup, "country");
    let dir = TempDir::new()?;
    let path = dir.path().join("mugration.nwk");

    write_tree_outputs(
      &AnnotatedTreeView::new(&graph)?,
      &btreemap! { TreeWriteKind::Nwk(NwkStyle::Beast) => path.clone() },
      CommandKind::Mugration,
      &NoopProgress,
    )?;

    let actual = read_file_to_string(&path)?;
    let expected = format!(
      "(A[&country={expected}]:0.1,B[&country=\"US\"]:0,C[&country=\"CH\"]:0.5)root[&country=\"US\"];\n"
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_mutations_precede_date_as_in_v0() -> Result<(), Report> {
    let setup = dated_setup()?;
    let mut edge_mutations = no_mutations(&setup.topology);
    edge_mutations.insert(
      parent_edge(&setup.topology, "A")?,
      vec![substitution(b'A', 54, b'G')?, substitution(b'T', 92, b'C')?],
    );
    let root_sequence = Seq::try_from_str("ACGT")?;
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
    let dir = TempDir::new()?;
    let path = dir.path().join("timetree.nwk");

    write_tree_outputs(
      &AnnotatedTreeView::new(&graph)?,
      &btreemap! { TreeWriteKind::Nwk(NwkStyle::Beast) => path.clone() },
      CommandKind::Timetree,
      &NoopProgress,
    )?;

    let actual = read_file_to_string(&path)?;
    let expected =
      "(A[&mutations=\"A55G,T93C\",date=2020.00]:0.1,B[&date=2021.00]:0,C[&date=2022.00]:0.5)root[&date=2023.00];\n";
    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::__tests__::test_tree_output::tests::helpers::{
      NUC_ALPHABET, annotated, no_mutations, node_key, parent_edge, topology_from,
    };
    use crate::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, TreeSequences};
    use crate::nwk_comments::nwk_node_comments;
    use eyre::Report;
    use treetime::seq::mutation::Mutation;
    use treetime_io::nwk::NewickValue;
    use treetime_primitives::Seq;

    pub(super) fn leaf_comments(mutations: Vec<Mutation>) -> Result<Vec<(String, NewickValue)>, Report> {
      let topology = topology_from("(A:0.1)root;")?;
      let mut edge_mutations = no_mutations(&topology);
      edge_mutations.insert(parent_edge(&topology, "A")?, mutations);
      let root_sequence = Seq::try_from_str("ACGT")?;
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
      let mut comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;
      Ok(comments.remove(&node_key(&topology, "A")).unwrap_or_default())
    }
  }
}
