use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
use crate::commands::prune::args::TreetimePruneArgs;
use crate::graph::edge::{GraphEdge, GraphEdgeKey, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey, Named};
use crate::gtr::get_gtr::{GtrModelName, JC69Params, get_gtr, jc69};
use crate::io::fasta::read_many_fasta;
use crate::io::nex::{NexWriteOptions, nex_write_file};
use crate::io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use crate::make_error;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use crate::utils::iter::iter_union;
use crate::utils::parse_delimited::{parse_delimited_file, parse_delimited_str};
use eyre::Report;
use itertools::Itertools;
use log::debug;
use maplit::{btreemap, btreeset};
use parking_lot::RwLock;
use serde::Serialize;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;

fn validate_args(args: &TreetimePruneArgs) -> Result<(), Report> {
  if args.prune_empty && args.input_fastas.is_empty() {
    return make_error!(
      "The --prune-empty requires --aln. Without sequence data, it's not possible to determine which branches lack mutations."
    );
  }

  Ok(())
}

pub fn run_prune(args: &TreetimePruneArgs) -> Result<(), Report> {
  validate_args(args)?;

  let TreetimePruneArgs {
    input_fastas,
    tree,
    alphabet,
    outdir,
    prune_short,
    prune_empty,
    prune_nodes_list,
    prune_nodes_list_delimiter,
    prune_nodes_list_file,
    prune_nodes_list_file_delimiter,
  } = args;

  let mut graph: GraphAncestral = nwk_read_file(tree)?;

  let partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = if *prune_empty {
    assert!(
      !input_fastas.is_empty(),
      "Input FASTA files are required for pruning empty branches."
    );
    let alphabet = Alphabet::new(alphabet.unwrap_or_default(), false)?;
    let aln = read_many_fasta(input_fastas, &alphabet)?;

    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
    // FIXME: spaghetti code: dummy gtr is replaced by real gtr here
    for partition in &partitions {
      let gtr = get_gtr(&GtrModelName::JC69, partition, &graph)?;
      partition.write_arc().gtr = gtr;
    }

    partitions
  } else {
    vec![]
  };

  let node_names = parse_node_names(
    prune_nodes_list.as_ref(),
    *prune_nodes_list_delimiter,
    prune_nodes_list_file.as_ref(),
    *prune_nodes_list_file_delimiter,
  )?;

  prune_nodes(&mut graph, &partitions, *prune_short, *prune_empty, &node_names)?;

  write_graph(outdir, &graph)?;

  Ok(())
}

fn parse_node_names(
  prune_nodes_list: Option<&String>,
  prune_nodes_list_delimiter: char,
  prune_nodes_list_file: Option<&PathBuf>,
  prune_nodes_list_file_delimiter: char,
) -> Result<BTreeSet<String>, Report> {
  let mut node_names = btreeset! {};

  if let Some(prune_nodes_list) = prune_nodes_list {
    node_names.extend(parse_delimited_str(prune_nodes_list, prune_nodes_list_delimiter as u8));
  }

  if let Some(prune_nodes_list_file) = prune_nodes_list_file {
    node_names.extend(parse_delimited_file(
      prune_nodes_list_file,
      prune_nodes_list_file_delimiter as u8,
    )?);
  }

  Ok(node_names)
}

fn prune_nodes(
  graph: &mut GraphAncestral,
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  prune_short: Option<f64>,
  prune_empty: bool,
  node_names: &BTreeSet<String>,
) -> Result<(), Report> {
  prune_internal_nodes(graph, partitions, prune_short, prune_empty, node_names)?;
  graph.build()?;
  prune_leaves(graph, partitions, node_names)?;
  graph.build()?;
  Ok(())
}

fn get_edge_num_muts(partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>, edge_key: GraphEdgeKey) -> Option<usize> {
  let mut total_muts = 0;
  let mut found_any = false;

  for partition in partitions {
    let partition = partition.read_arc();
    if let Some(edge_partition) = partition.edges.get(&edge_key) {
      total_muts += edge_partition.subs.len();
      found_any = true;
    }
  }

  found_any.then_some(total_muts)
}

fn prune_internal_nodes(
  graph: &mut GraphAncestral,
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  prune_short: Option<f64>,
  prune_empty: bool,
  node_names: &BTreeSet<String>,
) -> Result<(), Report> {
  #[allow(clippy::needless_collect)]
  let edges_to_collapse: Vec<_> = graph
    .get_edges()
    .iter()
    .filter_map(|edge| {
      let edge = edge.read_arc();
      let target_is_leaf = graph.is_leaf(edge.target());

      // Skip leaves in this pass
      if target_is_leaf {
        return None;
      }

      let weight = edge.payload().read_arc().weight();
      let should_prune_short = matches!((prune_short, weight), (Some(threshold), Some(weight)) if weight < threshold);

      let should_prune_empty = prune_empty && get_edge_num_muts(partitions, edge.key()) == Some(0);

      let target_node = graph.get_node(edge.target())?.read_arc().payload().read_arc();
      let name = target_node.name();
      let should_prune_by_name = name.is_some_and(|name| node_names.contains(name.as_ref()));

      let should_prune = should_prune_short || should_prune_empty || should_prune_by_name;
      should_prune.then(|| edge.key())
    })
    .collect();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing internal edge: {edge_key}");
    collapse_sparse_edge(graph, partitions, edge_key)
  })
}

fn prune_leaves(
  graph: &mut GraphAncestral,
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  node_names: &BTreeSet<String>,
) -> Result<(), Report> {
  #[allow(clippy::needless_collect)]
  let edges_to_collapse = graph
    .get_edges()
    .iter()
    .filter_map(|edge| {
      let edge = edge.read_arc();
      let target_is_leaf = graph.is_leaf(edge.target());

      if !target_is_leaf {
        return None;
      }

      let target_node = graph.get_node(edge.target())?.read_arc().payload().read_arc();
      let name = target_node.name();
      let should_prune_by_name = name.is_some_and(|name| node_names.contains(name.as_ref()));

      should_prune_by_name.then(|| edge.key())
    })
    .collect_vec();

  edges_to_collapse.into_iter().try_for_each(|edge_key| {
    debug!("Collapsing leaf edge: {edge_key}");
    collapse_sparse_edges_from_leaf_recursive(graph, partitions, edge_key)
  })
}

fn collapse_sparse_edges_from_leaf_recursive(
  graph: &mut GraphAncestral,
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  let mut current_edge_key = edge_key;

  loop {
    let parent_node_key = graph.get_source_node_key(current_edge_key)?;

    collapse_sparse_edge(graph, partitions, current_edge_key)?;

    let next_edge_key = if should_collapse_parent(graph, parent_node_key) {
      graph.parent_inbound_edge(parent_node_key)?
    } else {
      None
    };

    match next_edge_key {
      Some(key) => current_edge_key = key,
      None => break,
    }
  }

  Ok(())
}

fn should_collapse_parent(graph: &GraphAncestral, node_key: GraphNodeKey) -> bool {
  graph.has_at_most_one_child(node_key) && !graph.is_root(node_key)
}

fn collapse_sparse_edge(
  graph: &mut GraphAncestral,
  partitions: &Vec<Arc<RwLock<PartitionMarginalSparse>>>,
  edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  let (_, removed_edge, new_edges) = graph.collapse_edge(edge_key)?;
  let removed_edge = removed_edge.payload().read_arc();

  for new_edge in new_edges {
    let new_edge_key = new_edge.read_arc().key();
    let mut new_edge = new_edge.write_arc().payload().write_arc();

    // Sum branch lengths
    if let (Some(bl1), Some(bl2)) = (removed_edge.branch_length, new_edge.branch_length) {
      new_edge.branch_length = Some(bl1 + bl2);
    }

    // Union of substitutions per partition
    for partition in partitions {
      let mut partition = partition.write_arc();
      let removed_edge_subs = partition.edges[&edge_key].subs.clone(); // FIXME: avoid clone
      let new_edge = partition.edges.entry(new_edge_key).or_default();
      new_edge.subs = iter_union(&removed_edge_subs, &new_edge.subs).cloned().collect_vec();
    }
  }

  Ok(())
}

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + Serialize,
  E: GraphEdge + EdgeToNwk + Serialize,
  D: Send + Sync + Default + Serialize,
{
  nwk_write_file(
    outdir.as_ref().join("pruned_tree.nwk"),
    graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write_file(
    outdir.as_ref().join("pruned_tree.nexus"),
    graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::graph::graph::Graph;
  use crate::io::nwk::{nwk_read_str, nwk_write_str};
  use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
  use crate::representation::graph_sparse::SparseEdgePartition;
  use crate::seq::mutation::Sub;
  use pretty_assertions::assert_eq;

  fn create_test_graph_with_partitions(
    nwk: &str,
    edge_mutations: &[(usize, Option<usize>)], // (edge_index, num_mutations)
  ) -> Result<(GraphAncestral, Vec<Arc<RwLock<PartitionMarginalSparse>>>), Report> {
    let graph: GraphAncestral = nwk_read_str(nwk)?;

    let partitions = if edge_mutations.is_empty() {
      vec![]
    } else {
      let mut partition = PartitionMarginalSparse {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
        length: 100, // dummy length
        nodes: btreemap! {},
        edges: btreemap! {},
      };

      for (edge_index, num_muts) in edge_mutations {
        if let Some(edge) = graph.get_edges().get(*edge_index) {
          let edge_key = edge.read_arc().key();
          if let Some(num_muts) = num_muts {
            partition.edges.insert(
              edge_key,
              SparseEdgePartition {
                subs: (0..*num_muts).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
                ..Default::default()
              },
            );
          }
        }
      }

      vec![Arc::new(RwLock::new(partition))]
    };

    Ok((graph, partitions))
  }

  #[test]
  fn test_prune_nodes_basic() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.0,B:0.1)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.0), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(A:0[&mutations=\"\"],B:0.1[&mutations=\"\"])root[&mutations=\"\"];"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_with_threshold() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.01,B:0.02,C:0.1)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.05), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(A:0.01[&mutations=\"\"],B:0.02[&mutations=\"\"],C:0.1[&mutations=\"\"])root[&mutations=\"\"];"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_preserves_large_edges() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.1,B:0.2)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.01), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(A:0.1[&mutations=\"\"],B:0.2[&mutations=\"\"])root[&mutations=\"\"];"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_empty_graph() -> Result<(), Report> {
    let mut graph: GraphAncestral = Graph::new();
    let partitions = vec![];
    prune_nodes(&mut graph, &partitions, Some(0.0), false, &btreeset! {})?;
    assert!(graph.get_nodes().is_empty());
    Ok(())
  }

  #[test]
  fn test_prune_nodes_handles_none_weights() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.0,B:0.1)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.0), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(A:0[&mutations=\"\"],B:0.1[&mutations=\"\"])root[&mutations=\"\"];"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_preserves_terminal_nodes() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.00001,B:0.1)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.001), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(A:1.0e-5[&mutations=\"\"],B:0.1[&mutations=\"\"])root[&mutations=\"\"];"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_complex_tree() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions(
      "(((A:0,B:0.1)internal1:0.00002,(C:0.00003,D:0.1)internal2:0.1)internal3:0.00004,E:0.00005)root;",
      &[],
    )?;
    prune_nodes(&mut graph, &partitions, Some(0.01), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(E:5.0e-5[&mutations=\"\"],(C:3.0e-5[&mutations=\"\"],D:0.1[&mutations=\"\"])internal2:0.1[&mutations=\"\"],A:6.0e-5[&mutations=\"\"],B:0.1[&mutations=\"\"])root[&mutations=\"\"];"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_preserves_leaves() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      ..Default::default()
    });

    // Edge with no mutations to leaf A (should be preserved)
    graph.add_edge(
      root,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    // Edge with mutations to leaf B
    graph.add_edge(
      root,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_a_edge_key = graph.get_edges()[0].read_arc().key();
    let root_b_edge_key = graph.get_edges()[1].read_arc().key();

    // Edge to A has no mutations
    partition.edges.insert(root_a_edge_key, SparseEdgePartition::default());
    // Edge to B has 2 mutations
    partition.edges.insert(
      root_b_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..Default::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Both leaves should be preserved even if edge to A has no mutations
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_internal_nodes() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      ..Default::default()
    });

    // Internal edge with no mutations
    graph.add_edge(
      root,
      internal,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    // Leaf edges with mutations
    graph.add_edge(
      internal,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      internal,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_internal_edge_key = graph.get_edges()[0].read_arc().key();
    let internal_a_edge_key = graph.get_edges()[1].read_arc().key();
    let internal_b_edge_key = graph.get_edges()[2].read_arc().key();

    // Internal edge has no mutations
    partition
      .edges
      .insert(root_internal_edge_key, SparseEdgePartition::default());
    // Leaf edges have mutations
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..Default::default()
      },
    );
    partition.edges.insert(
      internal_b_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..Default::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Internal node should be collapsed, but leaves preserved
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_none_mutations() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });

    // Edge with None mutations (unknown number of mutations - should not be pruned)
    graph.add_edge(
      root,
      internal,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      internal,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let internal_a_edge_key = graph.get_edges()[1].read_arc().key();

    // Only add edge data for internal->A edge, leaving root->internal edge unknown (None mutations)
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..Default::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Internal node should be preserved when mutations is None (unknown)
    assert_eq!(graph.get_nodes().len(), 3); // root, internal, A
    assert_eq!(graph.get_edges().len(), 2); // root->internal, internal->A

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_simple_leaf_case() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });

    // Single leaf edge with no mutations should be preserved
    graph.add_edge(
      root,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_a_edge_key = graph.get_edges()[0].read_arc().key();
    partition.edges.insert(root_a_edge_key, SparseEdgePartition::default());

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // Leaf should be preserved even with no mutations
    assert_eq!(graph.get_nodes().len(), 2); // root and A
    assert_eq!(graph.get_edges().len(), 1); // root->A

    Ok(())
  }

  #[test]
  fn test_prune_nodes_combined_prune_short_and_empty() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let internal1 = graph.add_node(NodeAncestral {
      name: Some("internal1".to_owned()),
      ..Default::default()
    });
    let internal2 = graph.add_node(NodeAncestral {
      name: Some("internal2".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      ..Default::default()
    });

    // Short edge (should be pruned by threshold)
    graph.add_edge(
      root,
      internal1,
      EdgeAncestral {
        branch_length: Some(0.001),
        ..Default::default()
      },
    )?;
    // Empty edge (should be pruned by empty check)
    graph.add_edge(
      root,
      internal2,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    // Leaf edges (should be preserved)
    graph.add_edge(
      internal1,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      internal2,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_internal1_edge_key = graph.get_edges()[0].read_arc().key();
    let root_internal2_edge_key = graph.get_edges()[1].read_arc().key();
    let internal1_a_edge_key = graph.get_edges()[2].read_arc().key();
    let internal2_b_edge_key = graph.get_edges()[3].read_arc().key();

    // Short edge with mutations
    partition.edges.insert(
      root_internal1_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..Default::default()
      },
    );
    // Empty edge (no mutations)
    partition
      .edges
      .insert(root_internal2_edge_key, SparseEdgePartition::default());
    // Leaf edges
    partition.edges.insert(
      internal1_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..Default::default()
      },
    );
    partition.edges.insert(
      internal2_b_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..Default::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, Some(0.01), true, &btreeset! {})?;

    // Both internal nodes should be collapsed, leaves preserved
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2); // root->A, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_short_threshold_exact() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.05,B:0.05,C:0.051)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.05), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    assert_eq!(
      output_nwk,
      "(A:0.05[&mutations=\"\"],B:0.05[&mutations=\"\"],C:0.051[&mutations=\"\"])root[&mutations=\"\"];"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_short_threshold_below() -> Result<(), Report> {
    let (mut graph, partitions) = create_test_graph_with_partitions("(A:0.049,B:0.05,C:0.051)root;", &[])?;
    prune_nodes(&mut graph, &partitions, Some(0.05), false, &btreeset! {})?;
    let output_nwk = nwk_write_str(&graph, &NwkWriteOptions::default())?;
    // All edges remain because A, B, C are leaves and leaves are never collapsed
    assert_eq!(
      output_nwk,
      "(A:0.049[&mutations=\"\"],B:0.05[&mutations=\"\"],C:0.051[&mutations=\"\"])root[&mutations=\"\"];"
    );
    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_empty_complex_tree() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let internal1 = graph.add_node(NodeAncestral {
      name: Some("internal1".to_owned()),
      ..Default::default()
    });
    let internal2 = graph.add_node(NodeAncestral {
      name: Some("internal2".to_owned()),
      ..Default::default()
    });
    let internal3 = graph.add_node(NodeAncestral {
      name: Some("internal3".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      ..Default::default()
    });
    let c = graph.add_node(NodeAncestral {
      name: Some("C".to_owned()),
      ..Default::default()
    });
    let d = graph.add_node(NodeAncestral {
      name: Some("D".to_owned()),
      ..Default::default()
    });

    // Complex tree structure: root -> internal1 (with muts) -> A (leaf)
    //                              -> internal1 -> internal3 (no muts) -> C,D (leaves)
    //                         root -> internal2 (no muts) -> B (leaf)
    graph.add_edge(
      root,
      internal1,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?; // has muts
    graph.add_edge(
      root,
      internal2,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?; // no muts, to internal
    graph.add_edge(
      internal1,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?; // leaf, has muts
    graph.add_edge(
      internal1,
      internal3,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?; // no muts, to internal
    graph.add_edge(
      internal2,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?; // leaf, no muts (preserved)
    graph.add_edge(
      internal3,
      c,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?; // leaf, has muts
    graph.add_edge(
      internal3,
      d,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?; // leaf, has muts

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_internal1_edge_key = graph.get_edges()[0].read_arc().key();
    let root_internal2_edge_key = graph.get_edges()[1].read_arc().key();
    let internal1_a_edge_key = graph.get_edges()[2].read_arc().key();
    let internal1_internal3_edge_key = graph.get_edges()[3].read_arc().key();
    let internal2_b_edge_key = graph.get_edges()[4].read_arc().key();
    let internal3_c_edge_key = graph.get_edges()[5].read_arc().key();
    let internal3_d_edge_key = graph.get_edges()[6].read_arc().key();

    // Set up mutations for each edge
    partition.edges.insert(
      root_internal1_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..Default::default()
      },
    ); // has muts
    partition
      .edges
      .insert(root_internal2_edge_key, SparseEdgePartition::default()); // no muts, to internal
    partition.edges.insert(
      internal1_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..Default::default()
      },
    ); // leaf, has muts
    partition
      .edges
      .insert(internal1_internal3_edge_key, SparseEdgePartition::default()); // no muts, to internal
    partition
      .edges
      .insert(internal2_b_edge_key, SparseEdgePartition::default()); // leaf, no muts (preserved)
    partition.edges.insert(
      internal3_c_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..Default::default()
      },
    ); // leaf, has muts
    partition.edges.insert(
      internal3_d_edge_key,
      SparseEdgePartition {
        subs: (0_usize..2).map(|i| Sub::new('A', i, 'T').unwrap()).collect_vec(),
        ..Default::default()
      },
    ); // leaf, has muts

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, true, &btreeset! {})?;

    // internal2 and internal3 should be collapsed (empty internal edges), leaves preserved
    // Result: root -> internal1 -> A, C, D and root -> B
    assert_eq!(graph.get_nodes().len(), 6); // root, internal1, A, B, C, D
    assert_eq!(graph.get_edges().len(), 5); // root->internal1, internal1->A, internal1->C, internal1->D, root->B

    Ok(())
  }

  #[test]
  fn test_prune_nodes_prune_both_disabled() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let internal = graph.add_node(NodeAncestral {
      name: Some("internal".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });

    // Very short edge with no mutations - should be preserved when both pruning options disabled
    graph.add_edge(
      root,
      internal,
      EdgeAncestral {
        branch_length: Some(0.0001),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      internal,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;

    graph.build()?;

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(crate::alphabet::alphabet::AlphabetName::Nuc, false)?,
      length: 100,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    let root_internal_edge_key = graph.get_edges()[0].read_arc().key();
    let internal_a_edge_key = graph.get_edges()[1].read_arc().key();

    partition
      .edges
      .insert(root_internal_edge_key, SparseEdgePartition::default());
    partition.edges.insert(
      internal_a_edge_key,
      SparseEdgePartition {
        subs: vec![Sub::new('A', 0_usize, 'T').unwrap()],
        ..Default::default()
      },
    );

    let partitions = vec![Arc::new(RwLock::new(partition))];

    prune_nodes(&mut graph, &partitions, None, false, &btreeset! {})?;

    // Nothing should be collapsed
    assert_eq!(graph.get_nodes().len(), 3); // root, internal, A
    assert_eq!(graph.get_edges().len(), 2); // root->internal, internal->A

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_basic() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let internal1 = graph.add_node(NodeAncestral {
      name: Some("internal1".to_owned()),
      ..Default::default()
    });
    let internal2 = graph.add_node(NodeAncestral {
      name: Some("internal2".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      ..Default::default()
    });
    let c = graph.add_node(NodeAncestral {
      name: Some("C".to_owned()),
      ..Default::default()
    });

    // Tree structure: root -> internal1 -> internal2 -> A (the path to prune)
    //                     -> B (to keep)
    //                     -> C (to keep)
    graph.add_edge(
      root,
      internal1,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      root,
      b,
      EdgeAncestral {
        branch_length: Some(0.2),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      root,
      c,
      EdgeAncestral {
        branch_length: Some(0.3),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      internal1,
      internal2,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      internal2,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;

    graph.build()?;

    let partitions = vec![];

    // Find the edge leading to leaf A
    let a_inbound_edge = {
      let a_node = graph.get_node(a).unwrap();
      a_node.read_arc().inbound()[0]
    };

    // Recursively prune leaf A and its childless ancestors
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // The result should be: root -> B, root -> C (internal1 and internal2 should be removed)
    assert_eq!(graph.get_nodes().len(), 3); // root, B, C
    assert_eq!(graph.get_edges().len(), 2); // root->B, root->C

    // Verify the remaining nodes
    assert!(graph.get_node(root).is_some());
    assert!(graph.get_node(b).is_some());
    assert!(graph.get_node(c).is_some());

    // Verify removed nodes
    assert!(graph.get_node(a).is_none());
    assert!(graph.get_node(internal1).is_none());
    assert!(graph.get_node(internal2).is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_stops_at_node_with_children() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let internal1 = graph.add_node(NodeAncestral {
      name: Some("internal1".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });
    let b = graph.add_node(NodeAncestral {
      name: Some("B".to_owned()),
      ..Default::default()
    });

    // Tree structure: root -> internal1 -> A (to prune)
    //                               -> B (to keep)
    graph.add_edge(
      root,
      internal1,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      internal1,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.add_edge(
      internal1,
      b,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;

    graph.build()?;

    let partitions = vec![];

    // Find the edge leading to leaf A
    let a_inbound_edge = {
      let a_node = graph.get_node(a).unwrap();
      a_node.read_arc().inbound()[0]
    };

    // Recursively prune leaf A; internal1 becomes unary and should be collapsed upward
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // The result should be: root -> B (internal1 collapsed)
    assert_eq!(graph.get_nodes().len(), 2); // root, B
    assert_eq!(graph.get_edges().len(), 1); // root->B

    // Verify the remaining nodes
    assert!(graph.get_node(root).is_some());
    assert!(graph.get_node(b).is_some());

    // Verify removed node
    assert!(graph.get_node(a).is_none());
    assert!(graph.get_node(internal1).is_none());

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_stops_at_root() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });

    // Tree: root -> A (only child)
    graph.add_edge(
      root,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.build()?;

    let partitions = vec![];

    // Collapse the path starting at leaf A; should stop at root
    let a_inbound_edge = {
      let a_node = graph.get_node(a).unwrap();
      a_node.read_arc().inbound()[0]
    };
    collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, a_inbound_edge)?;

    // Only root should remain
    assert_eq!(graph.get_nodes().len(), 1);
    assert!(graph.get_node(root).is_some());
    assert!(graph.get_node(a).is_none());
    assert_eq!(graph.get_edges().len(), 0);

    Ok(())
  }

  #[test]
  fn test_collapse_sparse_edges_from_leaf_recursive_invalid_edge_key_errors() -> Result<(), Report> {
    let mut graph = GraphAncestral::new();

    let root = graph.add_node(NodeAncestral {
      name: Some("root".to_owned()),
      ..Default::default()
    });
    let a = graph.add_node(NodeAncestral {
      name: Some("A".to_owned()),
      ..Default::default()
    });

    graph.add_edge(
      root,
      a,
      EdgeAncestral {
        branch_length: Some(0.1),
        ..Default::default()
      },
    )?;
    graph.build()?;

    let partitions = vec![];

    // Use a non-existent edge key to ensure we surface an error
    let bogus = GraphEdgeKey(usize::MAX);
    let res = collapse_sparse_edges_from_leaf_recursive(&mut graph, &partitions, bogus);
    assert!(res.is_err());

    Ok(())
  }

  #[test]
  fn test_create_test_edge_num_muts_none_vs_some_zero() -> Result<(), Report> {
    // Test that we can distinguish between unknown mutations (None) and zero mutations (Some(0))
    let (graph, partitions) = create_test_graph_with_partitions("(A:0.1)root;", &[(0, None)])?;
    let edge_unknown_key = graph.get_edges()[0].read_arc().key();

    let (graph2, partitions2) = create_test_graph_with_partitions("(A:0.1)root;", &[(0, Some(0))])?;
    let edge_zero_key = graph2.get_edges()[0].read_arc().key();

    assert_eq!(get_edge_num_muts(&partitions, edge_unknown_key), None);
    assert_eq!(get_edge_num_muts(&partitions2, edge_zero_key), Some(0));

    Ok(())
  }
}
