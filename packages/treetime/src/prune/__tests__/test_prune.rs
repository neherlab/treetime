#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::Alphabet;

  use crate::prune::prune::prune_nodes;

  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs};

  use crate::seq::mutation::Sub;
  use crate::test_utils::find_edge_key;

  use eyre::Report;
  use itertools::Itertools;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::collections::BTreeSet;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
  use treetime_primitives::AsciiChar;
  use treetime_primitives::seq;

  use helpers::*;

  #[test]
  fn test_prune_nodes_basic() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.0,B:0.1)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.0),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_with_threshold() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.01,B:0.02,C:0.1)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.05),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.01,B:0.02,C:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_preserves_large_edges() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.1,B:0.2)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.01),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,B:0.2)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_empty_graph() -> Result<(), Report> {
    let mut graph = Graph::new();
    let names: BTreeMap<GraphNodeKey, Option<String>> = btreemap! {};
    let mut partitions = vec![];
    let mut branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>> = btreemap! {};
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.0),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    assert!(graph.get_nodes().next().is_none());
    Ok(())
  }

  #[test]
  fn test_prune_nodes_handles_none_weights() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.0,B:0.1)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.0),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_preserves_terminal_nodes() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.00001,B:0.1)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.001),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:1.0e-5,B:0.1)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_complex_tree() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) = create_test_graph_with_partitions(
      "(((A:0,B:0.1)internal1:0.00002,(C:0.00003,D:0.1)internal2:0.1)internal3:0.00004,E:0.00005)root;",
      &[],
    )?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.01),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(E:5.0e-5,(C:3.0e-5,D:0.1)internal2:0.1,A:6.00e-5,B:0.1)root;"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_preserves_leaves() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) = create_test_graph_with_named_edge_mutations(
      "(A:0.1,B:0.1)root;",
      &[("root", "A", None), ("root", "B", Some(2))],
    )?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      true,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_internal_nodes() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) = create_test_graph_with_named_edge_mutations(
      "((A:0.1,B:0.1)internal:0.1)root;",
      &[
        ("root", "internal", None),
        ("internal", "A", Some(1)),
        ("internal", "B", Some(2)),
      ],
    )?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      true,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_none_mutations() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_named_edge_mutations("((A:0.1)internal:0.1)root;", &[("internal", "A", Some(1))])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      true,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_simple_leaf_case() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_named_edge_mutations("(A:0.1)root;", &[("root", "A", None)])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      true,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 2);
    assert_eq!(graph.get_edges().count(), 1);

    Ok(())
  }

  #[test]
  fn test_prune_nodes_combined_prune_short_and_empty() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) = create_test_graph_with_named_edge_mutations(
      "((A:0.1)internal1:0.001,(B:0.1)internal2:0.1)root;",
      &[
        ("root", "internal1", Some(1)),
        ("root", "internal2", None),
        ("internal1", "A", None),
        ("internal2", "B", Some(2)),
      ],
    )?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.01),
      true,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_short_threshold_exact() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.05,B:0.05,C:0.051)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.05),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.05,B:0.05,C:0.051)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_short_threshold_below() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.049,B:0.05,C:0.051)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(0.05),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.049,B:0.05,C:0.051)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_complex_tree() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) = create_test_graph_with_named_edge_mutations(
      "(((C:0.1,D:0.1)internal3:0.1,A:0.1)internal1:0.1,(B:0.1)internal2:0.1)root;",
      &[
        ("root", "internal1", Some(2)),
        ("root", "internal2", Some(0)),
        ("internal1", "A", Some(1)),
        ("internal1", "internal3", Some(0)),
        ("internal2", "B", Some(0)),
        ("internal3", "C", Some(1)),
        ("internal3", "D", Some(2)),
      ],
    )?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      true,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 6);
    assert_eq!(graph.get_edges().count(), 5);

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_both_disabled() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) = create_test_graph_with_named_edge_mutations(
      "((A:0.1)internal:0.0001)root;",
      &[("root", "internal", Some(0)), ("internal", "A", Some(1))],
    )?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    Ok(())
  }

  #[test]
  fn test_prune_nodes_single_named_leaf() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.1,B:0.2,C:0.3)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "B".to_owned() },
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,C:0.3)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_multiple_named_leaves() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.1,B:0.2,C:0.3,D:0.4)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "A".to_owned(), "C".to_owned() },
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(B:0.2,D:0.4)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_nonexistent_name_is_noop() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.1,B:0.2)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "nonexistent".to_owned() },
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,B:0.2)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_named_internal_node() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("((A:0.1,B:0.2)internal:0.3,C:0.4)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(C:0.4,A:0.4,B:0.5)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_mixed_internal_and_leaf_names() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("((A:0.1,B:0.2)internal:0.3,C:0.4,D:0.5)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned(), "D".to_owned() },
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(C:0.4,A:0.4,B:0.5)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_empty_names_set_is_noop() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.1,B:0.2,C:0.3)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    let output_nwk = nwk_write_str(&graph, &names, &branch_lengths, &NwkWriteOptions::default())?;
    assert_eq!(output_nwk, "(A:0.1,B:0.2,C:0.3)root;");
    Ok(())
  }

  #[test]
  fn test_prune_nodes_all_leaves_preserves_root() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(A:0.1,B:0.2)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "A".to_owned(), "B".to_owned() },
      &names,
      &mut branch_lengths,
    )?;
    assert_eq!(graph.get_nodes().count(), 1);
    Ok(())
  }

  #[test]
  fn test_prune_nodes_deep_nested_leaf_removal() -> Result<(), Report> {
    let (mut graph, names, mut partitions, mut branch_lengths) =
      create_test_graph_with_partitions("(((A:0.1,B:0.2)i1:0.3,C:0.4)i2:0.5,D:0.6)root;", &[])?;
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "A".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 5);

    let leaf_names: BTreeSet<_> = graph
      .get_nodes()
      .filter(|n| graph.is_leaf(n.key()))
      .filter_map(|n| names.get(&n.key()).cloned().flatten())
      .collect();
    assert_eq!(leaf_names, btreeset! { "B".to_owned(), "C".to_owned(), "D".to_owned() });

    Ok(())
  }

  pub mod helpers {
    use super::*;

    pub fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub fn populate_test_nodes(partition: &mut PartitionMarginalSparse, graph: &Graph) {
      let alphabet = partition.alphabet.clone();
      let ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A')).take(partition.length).collect();
      if partition.root_sequence.is_empty() {
        partition.root_sequence = ref_seq;
      }
      for node in graph.get_nodes() {
        let key = node.key();
        partition
          .obs_nodes
          .entry(key)
          .or_insert_with(|| SparseNodeObs::empty(&alphabet));
      }
    }

    pub fn create_test_graph_with_partitions(
      nwk: &str,
      edge_mutations: &[(usize, Option<usize>)],
    ) -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        Vec<PartitionMarginalSparse>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(nwk)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;

      let partitions = if edge_mutations.is_empty() {
        vec![]
      } else {
        let mut partition = PartitionMarginalSparse {
          index: 0,
          alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
          length: 100,
          root_sequence: seq![],
          obs_nodes: btreemap! {},
          obs_edges: btreemap! {},
        };

        populate_test_nodes(&mut partition, &graph);

        for (edge_index, num_muts) in edge_mutations {
          if let Some(edge) = graph.get_edges().collect::<Vec<_>>().get(*edge_index) {
            let edge_key = edge.key();
            if let Some(num_muts) = num_muts {
              partition.obs_edges.insert(
                edge_key,
                SparseEdgeObs::with_fitch_subs(
                  (0..*num_muts)
                    .map(|i| Sub::new(c(b'A'), i, c(b'T')).unwrap())
                    .collect_vec(),
                ),
              );
            }
          }
        }

        vec![partition]
      };

      Ok((graph, names, partitions, branch_lengths))
    }

    pub fn create_test_graph_with_named_edge_mutations(
      nwk: &str,
      edge_mutations: &[(&str, &str, Option<usize>)],
    ) -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<String>>,
        Vec<PartitionMarginalSparse>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str(nwk)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;

      let partitions = if edge_mutations.is_empty() {
        vec![]
      } else {
        let mut partition = PartitionMarginalSparse {
          index: 0,
          alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
          length: 100,
          root_sequence: seq![],
          obs_nodes: btreemap! {},
          obs_edges: btreemap! {},
        };

        populate_test_nodes(&mut partition, &graph);

        for (source_name, target_name, num_muts) in edge_mutations {
          if let Some(edge_key) = find_edge_key(&graph, &names, source_name, target_name) {
            match num_muts {
              Some(n) => {
                partition.obs_edges.insert(
                  edge_key,
                  SparseEdgeObs::with_fitch_subs((0..*n).map(|i| Sub::new(c(b'A'), i, c(b'T')).unwrap()).collect_vec()),
                );
              },
              None => {
                partition.obs_edges.insert(edge_key, SparseEdgeObs::default());
              },
            }
          }
        }

        vec![partition]
      };

      Ok((graph, names, partitions, branch_lengths))
    }
  }
}
