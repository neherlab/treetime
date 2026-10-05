#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{c, edge_mutation_map};
  use app_output::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, Divergence, TreeDates, TreeSequences};
  use app_output::output_plan::{CommandKind, TreeWriteKind};
  use app_output::tree_output::write_tree_outputs;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use tempfile::TempDir;
  use treetime::progress::NoopProgress;
  use treetime::seq::mutation::Sub;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::{NwkStyle, nwk_read};
  use treetime_primitives::Seq;
  use treetime_utils::io::fs::read_file_to_string;

  #[test]
  fn test_timetree_nexus_output_includes_mutations_and_date() -> Result<(), Report> {
    let parsed = nwk_read(b"(A:0.1)root;".as_slice())?;
    let names = parsed.names();
    let graph = parsed.graph;
    let edge_subs = &[(
      0,
      vec![
        Sub::new(c(b'A'), 54_usize, c(b'G'))?,
        Sub::new(c(b'T'), 92_usize, c(b'C'))?,
      ],
    )];
    let edge_mutations = edge_mutation_map(&graph, edge_subs, &[]);
    let num_date: BTreeMap<GraphNodeKey, Option<f64>> = graph
      .get_nodes()
      .map(|node| (node.key(), node.is_leaf().then_some(2003.84)))
      .collect();
    let divergences: BTreeMap<GraphNodeKey, f64> = graph.get_nodes().map(|node| (node.key(), 0.0)).collect();
    let time_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = graph.get_edges().map(|edge| (edge.key(), None)).collect();
    let root_sequence = Seq::try_from_str("A")?;
    let excluded = BTreeSet::new();
    let annotated = AnnotatedGraph {
      graph: &graph,
      names: &names,
      divergence_branch_lengths: &time_lengths,
      time_branch_lengths: Some(&time_lengths),
      divergence: Divergence::Values(&divergences),
      branch_support: None,
      sequences: Some(TreeSequences {
        root_sequence: &root_sequence,
        edge_mutations: &edge_mutations,
        amino_acids: None,
      }),
      dates: Some(TreeDates {
        num_date: &num_date,
        confidence: None,
        excluded: &excluded,
      }),
      traits: None,
    };
    let dir = TempDir::new()?;
    let path = dir.path().join("timetree.nexus");

    write_tree_outputs(
      &AnnotatedTreeView::new(&annotated)?,
      &btreemap! { TreeWriteKind::Nexus(NwkStyle::Beast) => path.clone() },
      CommandKind::Timetree,
      &NoopProgress,
    )?;

    let actual = read_file_to_string(&path)?;
    let expected = indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=1;
        TaxLabels A;
      End;
      Begin Trees;
        Tree tree1=(A[&date=2003.84,mutations="A55G,T93C"])root;
      End;
    "#};
    assert_eq!(expected, actual);
    Ok(())
  }
}
