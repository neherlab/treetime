#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;

  use crate::o;
  use crate::optimize::topology::merge_shared_mutations::merge_shared_mutation_branches;
  use crate::prune::prune::{collapse_sparse_edges_from_leaf_recursive, get_edge_num_muts, prune_nodes};

  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs};
  use crate::pretty_assert_ulps_eq;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_relative_eq;
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
  use treetime_utils::make_report;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn populate_test_nodes(partition: &mut PartitionMarginalSparse, graph: &Graph) {
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

  fn create_test_graph_with_partitions(
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

  fn create_test_graph_with_named_edge_mutations(
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
  fn test_collapse_sparse_edges_from_leaf_recursive_basic() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((A:0.1)internal2:0.1)internal1:0.1,B:0.2,C:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mut partitions: Vec<PartitionMarginalSparse> = vec![];

    let a_inbound_edge =
      find_edge_key(&graph, &names, "internal2", "A").ok_or_else(|| make_report!("Edge internal2->A not found"))?;

    collapse_sparse_edges_from_leaf_recursive(&mut graph, &mut partitions, a_inbound_edge, &mut branch_lengths)?;

    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    assert!(find_node_key_by_name(&graph, &names, "root").is_some());
    assert!(find_node_key_by_name(&graph, &names, "B").is_some());
    assert!(find_node_key_by_name(&graph, &names, "C").is_some());

    assert!(find_node_key_by_name(&graph, &names, "A").is_none());
    assert!(find_node_key_by_name(&graph, &names, "internal1").is_none());
    assert!(find_node_key_by_name(&graph, &names, "internal2").is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_stops_at_node_with_children() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)internal1:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mut partitions: Vec<PartitionMarginalSparse> = vec![];

    let a_inbound_edge =
      find_edge_key(&graph, &names, "internal1", "A").ok_or_else(|| make_report!("Edge internal1->A not found"))?;

    collapse_sparse_edges_from_leaf_recursive(&mut graph, &mut partitions, a_inbound_edge, &mut branch_lengths)?;

    assert_eq!(graph.get_nodes().count(), 2);
    assert_eq!(graph.get_edges().count(), 1);

    assert!(find_node_key_by_name(&graph, &names, "root").is_some());
    assert!(find_node_key_by_name(&graph, &names, "B").is_some());

    assert!(find_node_key_by_name(&graph, &names, "A").is_none());
    assert!(find_node_key_by_name(&graph, &names, "internal1").is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_stops_at_root() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let mut partitions = vec![];

    let a_inbound_edge = find_edge_key(&graph, &names, "root", "A").unwrap();
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &mut partitions, a_inbound_edge, &mut branch_lengths)?;

    assert_eq!(graph.get_nodes().count(), 1);
    assert!(find_node_key_by_name(&graph, &names, "root").is_some());
    assert!(find_node_key_by_name(&graph, &names, "A").is_none());
    assert_eq!(graph.get_edges().count(), 0);

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_invalid_edge_key_errors() -> Result<(), Report> {
    let mut graph = Graph::new();

    let root = graph.add_node();
    let a = graph.add_node();

    let e_root_a = graph.add_edge(root, a)?;
    graph.build()?;

    let mut partitions = vec![];

    let bogus = GraphEdgeKey(usize::MAX);
    let mut branch_lengths = btreemap! { e_root_a => Some(0.1) };
    let res = collapse_sparse_edges_from_leaf_recursive(&mut graph, &mut partitions, bogus, &mut branch_lengths);
    assert!(res.is_err());

    Ok(())
  }

  #[test]
  fn test_create_test_edge_num_muts_none_vs_some_zero() -> Result<(), Report> {
    let (graph, names, partitions, _branch_lengths) = create_test_graph_with_partitions("(A:0.1)root;", &[(0, None)])?;
    let edge_unknown_key = graph.get_edges().collect::<Vec<_>>()[0].key();

    let (graph2, names2, partitions2, _branch_lengths2) =
      create_test_graph_with_partitions("(A:0.1)root;", &[(0, Some(0))])?;
    let edge_zero_key = graph2.get_edges().collect::<Vec<_>>()[0].key();

    assert_eq!(get_edge_num_muts(&partitions, edge_unknown_key)?, None);
    assert_eq!(get_edge_num_muts(&partitions2, edge_zero_key)?, Some(0));

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

  #[test]
  fn test_collapse_edge_compose_non_overlapping() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let root_internal_edge_key = find_edge_key(&graph, &names, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, &names, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, &names, "internal", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      root_sequence: seq![],
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
    };

    partition.obs_edges.insert(
      root_internal_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![
        Sub::new(c(b'A'), 0_usize, c(b'T'))?,
        Sub::new(c(b'A'), 1_usize, c(b'C'))?,
      ]),
    );
    partition.obs_edges.insert(
      internal_a_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![
        Sub::new(c(b'G'), 2_usize, c(b'T'))?,
        Sub::new(c(b'C'), 3_usize, c(b'A'))?,
      ]),
    );
    partition
      .obs_edges
      .insert(internal_b_edge_key, SparseEdgeObs::default());

    let mut partitions = vec![partition];

    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    let partition = &partitions[0];
    let mut a_muts_count = 0;
    let mut b_muts_count = 0;
    for edge in graph.get_edges() {
      let edge_key = edge.key();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());

      if let Some(edge_partition) = partition.obs_edges.get(&edge_key) {
        if target_name.as_deref() == Some("A") {
          a_muts_count = edge_partition.fitch_subs().len();
        } else if target_name.as_deref() == Some("B") {
          b_muts_count = edge_partition.fitch_subs().len();
        }
      }
    }

    assert_eq!(a_muts_count, 4);
    assert_eq!(b_muts_count, 2);

    Ok(())
  }

  #[test]
  fn test_collapse_edge_compose_chain() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let root_internal_edge_key = find_edge_key(&graph, &names, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, &names, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, &names, "internal", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      root_sequence: seq![],
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
    };

    partition.obs_edges.insert(
      root_internal_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 0_usize, c(b'G'))?]),
    );
    partition.obs_edges.insert(
      internal_a_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'G'), 0_usize, c(b'T'))?]),
    );
    partition
      .obs_edges
      .insert(internal_b_edge_key, SparseEdgeObs::default());

    let mut partitions = vec![partition];
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    let partition = &partitions[0];
    for edge in graph.get_edges() {
      let edge_key = edge.key();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());

      if target_name.as_deref() == Some("A") {
        let edge_partition = &partition.obs_edges[&edge_key];
        assert_eq!(edge_partition.fitch_subs().len(), 1);
        assert_eq!(edge_partition.fitch_subs()[0].reff(), c(b'A'));
        assert_eq!(edge_partition.fitch_subs()[0].qry(), c(b'T'));
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_compose_cancellation() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let root_internal_edge_key = find_edge_key(&graph, &names, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, &names, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, &names, "internal", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      root_sequence: seq![],
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
    };

    partition.obs_edges.insert(
      root_internal_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 0_usize, c(b'G'))?]),
    );
    partition.obs_edges.insert(
      internal_a_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'G'), 0_usize, c(b'A'))?]),
    );
    partition
      .obs_edges
      .insert(internal_b_edge_key, SparseEdgeObs::default());

    let mut partitions = vec![partition];
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    let partition = &partitions[0];
    for edge in graph.get_edges() {
      let edge_key = edge.key();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());

      if target_name.as_deref() == Some("A") {
        let edge_partition = &partition.obs_edges[&edge_key];
        assert_eq!(edge_partition.fitch_subs().len(), 0);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_compose_multiple_partitions() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let root_internal_edge_key = find_edge_key(&graph, &names, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, &names, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, &names, "internal", "B").unwrap();

    let mut partition1 = PartitionMarginalSparse {
      index: 0,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      root_sequence: seq![],
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
    };
    partition1.obs_edges.insert(
      root_internal_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?]),
    );
    partition1.obs_edges.insert(
      internal_a_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'G'), 1_usize, c(b'C'))?]),
    );
    partition1
      .obs_edges
      .insert(internal_b_edge_key, SparseEdgeObs::default());

    let mut partition2 = PartitionMarginalSparse {
      index: 1,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      root_sequence: seq![],
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
    };
    partition2.obs_edges.insert(
      root_internal_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![
        Sub::new(c(b'C'), 10_usize, c(b'A'))?,
        Sub::new(c(b'T'), 11_usize, c(b'G'))?,
      ]),
    );
    partition2.obs_edges.insert(
      internal_a_edge_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 12_usize, c(b'T'))?]),
    );
    partition2
      .obs_edges
      .insert(internal_b_edge_key, SparseEdgeObs::default());

    let mut partitions = vec![partition1, partition2];

    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    let p1 = &partitions[0];
    let p2 = &partitions[1];

    for edge in graph.get_edges() {
      let edge_key = edge.key();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());

      if target_name.as_deref() == Some("A") {
        assert_eq!(p1.obs_edges[&edge_key].fitch_subs().len(), 2);
        assert_eq!(p2.obs_edges[&edge_key].fitch_subs().len(), 3);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_sum_both_some() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.2,B:0.1)internal:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mut partitions = vec![];

    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    for edge in graph.get_edges() {
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());
      let branch_length = branch_lengths.get(&edge.key()).copied().flatten();

      if target_name.as_deref() == Some("A") {
        assert_relative_eq!(branch_length.unwrap(), 0.5, epsilon = 1e-6);
      } else if target_name.as_deref() == Some("B") {
        assert_relative_eq!(branch_length.unwrap(), 0.4, epsilon = 1e-6);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_sum_precision() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:2e-10,B:3e-10)internal:1e-10)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mut partitions = vec![];
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    for edge in graph.get_edges() {
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());
      let branch_length = branch_lengths.get(&edge.key()).copied().flatten();

      if target_name.as_deref() == Some("A") {
        pretty_assert_ulps_eq!(branch_length.unwrap(), 3e-10, max_ulps = 4);
      } else if target_name.as_deref() == Some("B") {
        pretty_assert_ulps_eq!(branch_length.unwrap(), 4e-10, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_none_plus_some() -> Result<(), Report> {
    let mut graph = Graph::new();

    let root = graph.add_node();
    let internal = graph.add_node();
    let a = graph.add_node();
    let b = graph.add_node();
    let names =
      btreemap! { root => Some(o!("root")), internal => Some(o!("internal")), a => Some(o!("A")), b => Some(o!("B")) };

    let e_root_internal = graph.add_edge(root, internal)?;
    let e_internal_a = graph.add_edge(internal, a)?;
    let e_internal_b = graph.add_edge(internal, b)?;

    graph.build()?;

    let mut partitions = vec![];
    let mut branch_lengths = btreemap! {
      e_root_internal => Some(0.5),
      e_internal_a => None,
      e_internal_b => Some(0.2),
    };
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    for edge in graph.get_edges() {
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());
      let branch_length = branch_lengths.get(&edge.key()).copied().flatten();

      if target_name.as_deref() == Some("A") {
        assert!(branch_length.is_none());
      } else if target_name.as_deref() == Some("B") {
        pretty_assert_ulps_eq!(branch_length.unwrap(), 0.7, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_some_plus_none() -> Result<(), Report> {
    let mut graph = Graph::new();

    let root = graph.add_node();
    let internal = graph.add_node();
    let a = graph.add_node();
    let b = graph.add_node();
    let names =
      btreemap! { root => Some(o!("root")), internal => Some(o!("internal")), a => Some(o!("A")), b => Some(o!("B")) };

    let e_root_internal = graph.add_edge(root, internal)?;
    let e_internal_a = graph.add_edge(internal, a)?;
    let e_internal_b = graph.add_edge(internal, b)?;

    graph.build()?;

    let mut partitions = vec![];
    let mut branch_lengths = btreemap! {
      e_root_internal => None,
      e_internal_a => Some(0.3),
      e_internal_b => Some(0.2),
    };
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    for edge in graph.get_edges() {
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());
      let branch_length = branch_lengths.get(&edge.key()).copied().flatten();

      if target_name.as_deref() == Some("A") {
        pretty_assert_ulps_eq!(branch_length.unwrap(), 0.3, max_ulps = 4);
      } else if target_name.as_deref() == Some("B") {
        pretty_assert_ulps_eq!(branch_length.unwrap(), 0.2, max_ulps = 4);
      }
    }

    Ok(())
  }

  #[test]
  fn test_collapse_edge_branch_length_both_none() -> Result<(), Report> {
    let mut graph = Graph::new();

    let root = graph.add_node();
    let internal = graph.add_node();
    let a = graph.add_node();
    let b = graph.add_node();
    let names =
      btreemap! { root => Some(o!("root")), internal => Some(o!("internal")), a => Some(o!("A")), b => Some(o!("B")) };

    let e_root_internal = graph.add_edge(root, internal)?;
    let e_internal_a = graph.add_edge(internal, a)?;
    let e_internal_b = graph.add_edge(internal, b)?;

    graph.build()?;

    let mut partitions = vec![];
    let mut branch_lengths = btreemap! {
      e_root_internal => None,
      e_internal_a => None,
      e_internal_b => None,
    };
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    for edge in graph.get_edges() {
      let branch_length = branch_lengths.get(&edge.key()).copied().flatten();
      assert!(branch_length.is_none());
    }

    Ok(())
  }

  #[test]
  fn test_prune_then_merge_exposes_hidden_polytomy() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:1e-8,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      root_sequence: seq![],
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
    };

    populate_test_nodes(&mut partition, &graph);

    let pi_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    partition.obs_edges.insert(pi_key, SparseEdgeObs::default());

    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    partition.obs_edges.insert(
      ia_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?]),
    );

    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();
    partition.obs_edges.insert(
      ib_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?]),
    );

    let pc_key = find_edge_key(&graph, &names, "root", "C").unwrap();
    partition.obs_edges.insert(
      pc_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 0_usize, c(b'T'))?]),
    );

    let pd_key = find_edge_key(&graph, &names, "root", "D").unwrap();
    partition.obs_edges.insert(
      pd_key,
      SparseEdgeObs::with_fitch_subs(vec![Sub::new(c(b'A'), 5_usize, c(b'T'))?]),
    );

    let mut partitions = vec![partition];

    prune_nodes(
      &mut graph,
      &mut partitions,
      Some(1e-6),
      false,
      &btreeset! {},
      &names,
      &mut branch_lengths,
    )?;
    assert!(
      find_node_key_by_name(&graph, &names, "I").is_none(),
      "I should be collapsed by prune"
    );

    let merged = merge_shared_mutation_branches(&mut graph, &mut partitions, &mut branch_lengths)?;
    graph.build()?;

    assert_eq!(merged, 1, "one group merge should place A, B, C under one node");

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.degree_out(), 2);

    assert!(find_node_key_by_name(&graph, &names, "D").is_some());

    Ok(())
  }

  #[test]
  fn test_collapse_edge_indel_preservation() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)internal:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let root_internal_edge_key = find_edge_key(&graph, &names, "root", "internal").unwrap();
    let internal_a_edge_key = find_edge_key(&graph, &names, "internal", "A").unwrap();
    let internal_b_edge_key = find_edge_key(&graph, &names, "internal", "B").unwrap();

    let mut partition = PartitionMarginalSparse {
      index: 0,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc)?,
      length: 100,
      root_sequence: seq![],
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
    };

    let parent_indel = InDel::del((10, 15), [c(b'A'), c(b'C'), c(b'G'), c(b'T'), c(b'A')].as_slice())?;
    let child_indel = InDel::ins((20, 23), [c(b'G'), c(b'G'), c(b'C')].as_slice())?;

    partition.obs_edges.insert(
      root_internal_edge_key,
      SparseEdgeObs::with_fitch_subs_and_indels(vec![], vec![parent_indel]),
    );
    partition.obs_edges.insert(
      internal_a_edge_key,
      SparseEdgeObs::with_fitch_subs_and_indels(vec![], vec![child_indel]),
    );
    partition
      .obs_edges
      .insert(internal_b_edge_key, SparseEdgeObs::default());

    let mut partitions = vec![partition];
    prune_nodes(
      &mut graph,
      &mut partitions,
      None,
      false,
      &btreeset! { "internal".to_owned() },
      &names,
      &mut branch_lengths,
    )?;

    let partition = &partitions[0];
    for edge in graph.get_edges() {
      let edge_key = edge.key();
      let target_key = edge.target();
      let target_name = graph
        .get_node(target_key)
        .and_then(|n| names.get(&n.key()).cloned().flatten());

      if target_name.as_deref() == Some("A") {
        let edge_partition = &partition.obs_edges[&edge_key];
        assert_eq!(edge_partition.indels.len(), 2);
        assert_eq!(edge_partition.indels[0].range, (10, 15));
        assert!(edge_partition.indels[0].is_deletion());
        assert_eq!(edge_partition.indels[1].range, (20, 23));
        assert!(!edge_partition.indels[1].is_deletion());
      } else if target_name.as_deref() == Some("B") {
        let edge_partition = &partition.obs_edges[&edge_key];
        assert_eq!(edge_partition.indels.len(), 1);
        assert_eq!(edge_partition.indels[0].range, (10, 15));
      }
    }

    Ok(())
  }
}
