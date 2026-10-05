#[cfg(test)]
mod tests {
  use crate::__tests__::test_tree_output::tests::helpers::{
    annotated, dated_graph, dated_setup, mugration_graph, mugration_setup, no_mutations, node_key, parent_edge,
    substitution, topology_from,
  };
  use crate::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, TreeSequences};
  use crate::nwk_comments::nwk_node_comments;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack};
  use treetime_primitives::Seq;
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

    assert_eq!(btreemap! { o!("mutations") => o!("A1T,C2-,G3-,G6C") }, actual);
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

    assert_eq!(btreemap! { o!("mutations") => o!("A11T,G31C,C51G") }, actual);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_edge_without_mutations_has_none() -> Result<(), Report> {
    assert_eq!(BTreeMap::new(), helpers::leaf_comments(vec![])?);
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
        root_sequence: &root_sequence,
        edge_mutations: &edge_mutations,
        amino_acids: None,
      }),
      ..annotated(&topology)
    };

    let comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;

    assert_eq!(BTreeMap::new(), comments[&node_key(&topology, "root")]);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_dates_without_sequences_give_date_comments() -> Result<(), Report> {
    let setup = dated_setup()?;
    let graph = dated_graph(&setup, None);

    let comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;

    let actual: BTreeMap<String, BTreeMap<String, String>> = comments
      .into_iter()
      .map(|(key, comments)| (setup.topology.names[&key].clone().unwrap(), comments))
      .collect();
    let expected = btreemap! {
      o!("A") => btreemap! { o!("date") => o!("2020.00") },
      o!("B") => btreemap! { o!("date") => o!("2021.00") },
      o!("C") => btreemap! { o!("date") => o!("2022.00") },
      o!("root") => btreemap! { o!("date") => o!("2023.00") },
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
        root_sequence: &root_sequence,
        edge_mutations: &edge_mutations,
        amino_acids: None,
      }),
      ..dated_graph(&setup, None)
    };

    let comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;

    let expected = btreemap! { o!("date") => o!("2020.00"), o!("mutations") => o!("A1T") };
    assert_eq!(expected, comments[&node_key(&setup.topology, "A")]);
    Ok(())
  }

  #[test]
  fn test_nwk_comments_traits_give_trait_comments() -> Result<(), Report> {
    let setup = mugration_setup()?;
    let support = btreemap! {};
    let graph = mugration_graph(&setup, "country", &support);

    let comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;

    let actual: BTreeMap<String, BTreeMap<String, String>> = comments
      .into_iter()
      .map(|(key, comments)| (setup.topology.names[&key].clone().unwrap(), comments))
      .collect();
    let expected = btreemap! {
      o!("A") => btreemap! { o!("country") => o!("CH") },
      o!("B") => btreemap! { o!("country") => o!("US") },
      o!("C") => btreemap! { o!("country") => o!("CH") },
      o!("root") => btreemap! { o!("country") => o!("US") },
    };
    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::__tests__::test_tree_output::tests::helpers::{
      annotated, no_mutations, node_key, parent_edge, topology_from,
    };
    use crate::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, TreeSequences};
    use crate::nwk_comments::nwk_node_comments;
    use eyre::Report;
    use std::collections::BTreeMap;
    use treetime::seq::mutation::Mutation;
    use treetime_primitives::Seq;

    pub(super) fn leaf_comments(mutations: Vec<Mutation>) -> Result<BTreeMap<String, String>, Report> {
      let topology = topology_from("(A:0.1)root;")?;
      let mut edge_mutations = no_mutations(&topology);
      edge_mutations.insert(parent_edge(&topology, "A")?, mutations);
      let root_sequence = Seq::try_from_str("ACGT")?;
      let graph = AnnotatedGraph {
        sequences: Some(TreeSequences {
          root_sequence: &root_sequence,
          edge_mutations: &edge_mutations,
          amino_acids: None,
        }),
        ..annotated(&topology)
      };
      let mut comments = nwk_node_comments(&AnnotatedTreeView::new(&graph)?)?;
      Ok(comments.remove(&node_key(&topology, "A")).unwrap_or_default())
    }
  }
}
