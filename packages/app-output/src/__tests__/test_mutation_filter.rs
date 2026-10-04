#[cfg(test)]
mod tests {
  use crate::mutation_filter::UnknownMutationFilter;
  use eyre::Report;
  use helpers::{edge_strings, node_events, nucleotide_edges, tree};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime::ancestral::aa::AaNodeData;
  use treetime::seq::mutation::{AlignedMutation, MutationEvent};
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::{o, vec_of_owned};

  #[rustfmt::skip]
  #[rstest]
  #[case::default_drops_unknown( false, vec_of_owned!["C3K", "A4T"])]
  #[case::report_ambiguous(      true,  vec_of_owned!["C1N", "C3K", "A4T"])]
  #[trace]
  fn test_mutation_filter_leaf_unknown(#[case] report_ambiguous: bool, #[case] expected_b: Vec<String>) -> Result<(), Report> {
    let (graph, names) = tree()?;
    let edges = nucleotide_edges(&graph, &names, &btreemap! { "B" => vec!["C1N", "C3K", "A4T"] })?;

    let actual = UnknownMutationFilter::new(AsciiChar::try_new(b'N')?, report_ambiguous).reported_edge_mutations(&graph, edges)?;

    let expected = btreemap! { o!("B") => expected_b };
    assert_eq!(expected, edge_strings(&graph, &names, &actual)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default_bridges_unknown_node(false, btreemap! {
    o!("X") => vec![],               o!("A") => vec_of_owned!["A2G"], o!("B") => vec![],
  })]
  #[case::report_ambiguous(            true,  btreemap! {
    o!("X") => vec_of_owned!["A2N"], o!("A") => vec_of_owned!["N2G"], o!("B") => vec_of_owned!["N2A"],
  })]
  #[trace]
  fn test_mutation_filter_internal_unknown(
    #[case] report_ambiguous: bool,
    #[case] expected: BTreeMap<String, Vec<String>>,
  ) -> Result<(), Report> {
    let (graph, names) = tree()?;
    let edges = nucleotide_edges(&graph, &names, &btreemap! {
      "X" => vec!["A2N"],
      "A" => vec!["N2G"],
      "B" => vec!["N2A"],
    })?;

    let actual = UnknownMutationFilter::new(AsciiChar::try_new(b'N')?, report_ambiguous).reported_edge_mutations(&graph, edges)?;

    assert_eq!(expected, edge_strings(&graph, &names, &actual)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::leaf_edges(        btreemap! { "A" => vec!["C1N", "G2T"], "B" => vec!["C1N"] },
                             btreemap! { "A" => vec!["G2T"],        "B" => vec![] })]
  #[case::below_unknown_node(btreemap! { "X" => vec!["A2N"], "A" => vec!["N2G", "C3N"], "B" => vec!["N2A", "T4N"] },
                             btreemap! { "X" => vec!["A2N"], "A" => vec!["N2G"],        "B" => vec!["N2A"] })]
  #[trace]
  fn test_mutation_filter_hidden_unknown_ignores_leaf_subs_into_unknown(
    #[case] with_leaf_unknown: BTreeMap<&str, Vec<&str>>,
    #[case] without_leaf_unknown: BTreeMap<&str, Vec<&str>>,
  ) -> Result<(), Report> {
    let (graph, names) = tree()?;
    let filter = UnknownMutationFilter::hiding_unknown(AsciiChar::try_new(b'N')?);

    let expected = filter.reported_edge_mutations(&graph, nucleotide_edges(&graph, &names, &with_leaf_unknown)?)?;
    let actual = filter.reported_edge_mutations(&graph, nucleotide_edges(&graph, &names, &without_leaf_unknown)?)?;

    assert_eq!(edge_strings(&graph, &names, &expected)?, edge_strings(&graph, &names, &actual)?);
    Ok(())
  }

  #[test]
  fn test_mutation_filter_keeps_indels() -> Result<(), Report> {
    let (graph, names) = tree()?;
    let deletion = MutationEvent::Deletion(AlignedMutation {
      range: (4, 6),
      sequence: Seq::try_from_str("AC")?,
    });
    let mut edges = nucleotide_edges(&graph, &names, &btreemap! { "C" => vec!["C1N"] })?;
    let c_edge = helpers::edge_key(&graph, &names, "C")?;
    edges.entry(c_edge).or_default().push(helpers::nucleotide(deletion));

    let actual =
      UnknownMutationFilter::hiding_unknown(AsciiChar::try_new(b'N')?).reported_edge_mutations(&graph, edges)?;

    let expected = btreemap! { o!("C") => vec_of_owned!["A5-", "C6-"] };
    assert_eq!(expected, edge_strings(&graph, &names, &actual)?);
    Ok(())
  }

  #[test]
  fn test_mutation_filter_amino_acid_bridges_unknown_per_cds() -> Result<(), Report> {
    let (graph, names) = tree()?;
    let node_data = AaNodeData {
      node_aa_mutations: node_events(
        &names,
        &btreemap! {
          "root" => btreemap! { "S" => vec!["K7R"] },
          "X" => btreemap! { "S" => vec!["A1X"], "E" => vec!["M3X"] },
          "A" => btreemap! { "S" => vec!["X1K"], "E" => vec!["X3M"] },
          "B" => btreemap! { "S" => vec!["X1A", "K2B"], "E" => vec![] },
          "C" => btreemap! { "S" => vec!["A1X"], "E" => vec![] },
        },
      )?,
      ..AaNodeData::default()
    };

    let actual =
      UnknownMutationFilter::hiding_unknown(AsciiChar::try_new(b'X')?).reported_aa_node_data(&graph, node_data)?;

    let expected = node_events(
      &names,
      &btreemap! {
        "root" => btreemap! { "S" => vec!["K7R"] },
        "X" => btreemap! { "S" => vec![], "E" => vec![] },
        "A" => btreemap! { "S" => vec!["A1K"], "E" => vec![] },
        "B" => btreemap! { "S" => vec!["K2B"], "E" => vec![] },
        "C" => btreemap! { "S" => vec![], "E" => vec![] },
      },
    )?;
    assert_eq!(expected, actual.node_aa_mutations);
    Ok(())
  }

  mod helpers {
    use eyre::Report;
    use std::collections::BTreeMap;
    use std::str::FromStr;
    use treetime::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub, mutation_event_strings};
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read_str;

    pub(super) type Names = BTreeMap<GraphNodeKey, Option<String>>;

    pub(super) fn tree() -> Result<(Graph, Names), Report> {
      let parsed = nwk_read_str("((A:1,B:1)X:1,C:1)root;")?;
      let names = parsed.names();
      Ok((parsed.graph, names))
    }

    pub(super) fn node_key(names: &Names, name: &str) -> GraphNodeKey {
      names
        .iter()
        .find(|(_, node_name)| node_name.as_deref() == Some(name))
        .map(|(&key, _)| key)
        .expect("node must exist")
    }

    pub(super) fn edge_key(graph: &Graph, names: &Names, name: &str) -> Result<GraphEdgeKey, Report> {
      let (_, edge_key) = graph
        .node_parent(node_key(names, name))?
        .expect("node must have a parent");
      Ok(edge_key)
    }

    pub(super) fn nucleotide(event: MutationEvent) -> Mutation {
      Mutation {
        track: MutationTrack::Nucleotide,
        event,
      }
    }

    pub(super) fn events(strings: &[&str]) -> Result<Vec<MutationEvent>, Report> {
      strings
        .iter()
        .map(|string| Ok(MutationEvent::Substitution(Sub::from_str(string)?)))
        .collect()
    }

    pub(super) fn nucleotide_edges(
      graph: &Graph,
      names: &Names,
      mutations: &BTreeMap<&str, Vec<&str>>,
    ) -> Result<BTreeMap<GraphEdgeKey, Vec<Mutation>>, Report> {
      mutations
        .iter()
        .map(|(name, strings)| {
          let mutations = events(strings)?.into_iter().map(nucleotide).collect();
          Ok((edge_key(graph, names, name)?, mutations))
        })
        .collect()
    }

    pub(super) fn node_events(
      names: &Names,
      mutations: &BTreeMap<&str, BTreeMap<&str, Vec<&str>>>,
    ) -> Result<BTreeMap<GraphNodeKey, BTreeMap<String, Vec<MutationEvent>>>, Report> {
      mutations
        .iter()
        .map(|(name, cds_strings)| {
          let cds_events = cds_strings
            .iter()
            .map(|(cds, strings)| Ok(((*cds).to_owned(), events(strings)?)))
            .collect::<Result<_, Report>>()?;
          Ok((node_key(names, name), cds_events))
        })
        .collect()
    }

    pub(super) fn edge_strings(
      graph: &Graph,
      names: &Names,
      edges: &BTreeMap<GraphEdgeKey, Vec<Mutation>>,
    ) -> Result<BTreeMap<String, Vec<String>>, Report> {
      edges
        .iter()
        .map(|(&edge_key, mutations)| {
          let child = graph.get_edge(edge_key).expect("edge must exist").target();
          let strings = mutations
            .iter()
            .map(|mutation| mutation_event_strings(&mutation.event))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .flatten()
            .collect();
          Ok((names[&child].clone().expect("every node is named"), strings))
        })
        .collect()
    }
  }
}
