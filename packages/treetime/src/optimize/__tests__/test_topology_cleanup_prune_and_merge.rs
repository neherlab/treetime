#[cfg(test)]
mod tests {

  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};

  use crate::optimize::params::TopologyOps;
  use crate::optimize::run_loop::prune_and_merge_in_loop;

  use crate::partition::storage::sparse::SparseEdgeObs;

  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use eyre::Report;

  use pretty_assertions::assert_eq;

  use treetime_graph::graph::Graph;

  use treetime_io::nwk::nwk_read_str;

  use crate::optimize::__tests__::test_topology_cleanup::tests::helpers::*;

  #[test]
  fn test_optimize_prune_and_merge_empty_list() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let sparse: Vec<SparseReconstruction> = vec![];
    let dense: Vec<DenseReconstruction> = vec![];

    let mut names_tt_13 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_13,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(!changed);
    assert_eq!(graph.get_nodes().count(), 4);
    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_collapses_and_merges() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut partition = empty_sparse_recon()?;

    populate_test_nodes(&mut partition, &graph);

    partition.partition.obs_edges.insert(ri_key, SparseEdgeObs::default());

    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();
    let rc_key = find_edge_key(&graph, &names, "root", "C").unwrap();
    let rd_key = find_edge_key(&graph, &names, "root", "D").unwrap();

    partition
      .partition
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(ib_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rc_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rd_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    branch_lengths.insert(ri_key, Some(0.0));

    let mut names_tt_12 = names.clone();
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[ri_key],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_12,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed);

    assert!(find_node_key_by_name(&graph, &names, "I").is_none());

    assert!(find_node_key_by_name(&graph, &names, "D").is_some());

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.degree_out(), 2);

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_hoists_reversion_without_collapse() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mut partition = empty_sparse_recon()?;
    populate_test_nodes(&mut partition, &graph);

    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let vc1 = find_edge_key(&graph, &names, "V", "C1").unwrap();
    let vc2 = find_edge_key(&graph, &names, "V", "C2").unwrap();
    let vc3 = find_edge_key(&graph, &names, "V", "C3").unwrap();
    partition.partition.obs_edges.insert(
      uv,
      SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
    );
    partition
      .partition
      .obs_edges
      .insert(vc1, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition
      .partition
      .obs_edges
      .insert(vc2, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition.partition.obs_edges.insert(vc3, SparseEdgeObs::default());

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    let mut names_tt_9 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_9,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed, "reversion polytomy must be resolved even without a collapse");

    let p = &sparse[0];
    let total_subs: usize = graph
      .get_edges()
      .filter_map(|e| p.partition.obs_edges.get(&e.key()))
      .map(|e| e.fitch_subs().len())
      .sum();
    assert_eq!(total_subs, 2, "reaches the parsimony optimum");

    let reversion_remains = graph
      .get_edges()
      .filter_map(|e| p.partition.obs_edges.get(&e.key()))
      .any(|e| e.fitch_subs().contains(&sub(b'T', 0, b'A')));
    assert!(!reversion_remains, "reversion must be removed");

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_names_new_nodes() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut partition = empty_sparse_recon()?;

    populate_test_nodes(&mut partition, &graph);

    partition.partition.obs_edges.insert(ri_key, SparseEdgeObs::default());

    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();
    let rc_key = find_edge_key(&graph, &names, "root", "C").unwrap();
    let rd_key = find_edge_key(&graph, &names, "root", "D").unwrap();

    partition
      .partition
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(ib_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rc_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rd_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    branch_lengths.insert(ri_key, Some(0.0));

    let mut names_tt_6 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[ri_key],
      TopologyOps::default(),
      &mut branch_lengths,
      &mut names_tt_6,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed);

    let mut names: Vec<String> = graph
      .get_nodes()
      .filter_map(|n| names_tt_6.get(&n.key()).cloned().flatten())
      .collect();
    names.sort();

    assert_eq!(names, vec!["A", "B", "C", "D", "NODE_0000000", "root"]);

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_merge_disabled_keeps_polytomy() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0,C:0.1,D:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut partition = empty_sparse_recon()?;
    populate_test_nodes(&mut partition, &graph);
    partition.partition.obs_edges.insert(ri_key, SparseEdgeObs::default());

    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();
    let rc_key = find_edge_key(&graph, &names, "root", "C").unwrap();
    let rd_key = find_edge_key(&graph, &names, "root", "D").unwrap();
    partition
      .partition
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(ib_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rc_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .partition
      .obs_edges
      .insert(rd_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    branch_lengths.insert(ri_key, Some(0.0));

    let ops = TopologyOps {
      merge_siblings: false,
      ..TopologyOps::default()
    };
    let mut names_tt_5 = names.clone();
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[ri_key],
      ops,
      &mut branch_lengths,
      &mut names_tt_5,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed, "collapse still fires even with merge disabled");

    assert!(find_node_key_by_name(&graph, &names, "I").is_none());

    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_node = graph.get_node(root_key).unwrap();
    assert_eq!(root_node.degree_out(), 4);

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_flip_disabled_keeps_reversion() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mut partition = empty_sparse_recon()?;
    populate_test_nodes(&mut partition, &graph);

    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let vc1 = find_edge_key(&graph, &names, "V", "C1").unwrap();
    let vc2 = find_edge_key(&graph, &names, "V", "C2").unwrap();
    let vc3 = find_edge_key(&graph, &names, "V", "C3").unwrap();
    partition.partition.obs_edges.insert(
      uv,
      SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
    );
    partition
      .partition
      .obs_edges
      .insert(vc1, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition
      .partition
      .obs_edges
      .insert(vc2, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition.partition.obs_edges.insert(vc3, SparseEdgeObs::default());

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    let ops = TopologyOps {
      flip_parent_child: false,
      ..TopologyOps::default()
    };
    let mut names_tt_4 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[],
      ops,
      &mut branch_lengths,
      &mut names_tt_4,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(changed, "merge still groups the reverting siblings");

    let p = &sparse[0];
    let reversion_remains = graph
      .get_edges()
      .filter_map(|e| p.partition.obs_edges.get(&e.key()))
      .any(|e| e.fitch_subs().contains(&sub(b'T', 0, b'A')));
    assert!(
      reversion_remains,
      "reversion is kept when flip-parent-child is disabled"
    );

    Ok(())
  }

  #[test]
  fn test_optimize_prune_and_merge_all_ops_disabled_is_noop() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(((C1:0.1,C2:0.1,C3:0.1)V:0.2)U:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let mut partition = empty_sparse_recon()?;
    populate_test_nodes(&mut partition, &graph);

    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let vc1 = find_edge_key(&graph, &names, "V", "C1").unwrap();
    let vc2 = find_edge_key(&graph, &names, "V", "C2").unwrap();
    let vc3 = find_edge_key(&graph, &names, "V", "C3").unwrap();
    partition.partition.obs_edges.insert(
      uv,
      SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G')]),
    );
    partition
      .partition
      .obs_edges
      .insert(vc1, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition
      .partition
      .obs_edges
      .insert(vc2, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));
    partition.partition.obs_edges.insert(vc3, SparseEdgeObs::default());

    let sparse = vec![partition];
    let dense: Vec<DenseReconstruction> = vec![];

    let node_count_before = graph.get_nodes().count();
    let ops = TopologyOps {
      collapse_short_branches: false,
      merge_siblings: false,
      flip_parent_child: false,
    };
    let mut names_tt_3 = names;
    let cleanup = prune_and_merge_in_loop(
      &mut graph,
      sparse,
      dense,
      &[],
      ops,
      &mut branch_lengths,
      &mut names_tt_3,
    )?;
    let sparse = cleanup.sparse_partitions;
    let dense = cleanup.dense_partitions;
    let changed = cleanup.topology_changed;
    assert!(!changed, "no topology step runs when all are disabled");
    assert_eq!(graph.get_nodes().count(), node_count_before);

    let p = &sparse[0];
    let total_subs: usize = graph
      .get_edges()
      .filter_map(|e| p.partition.obs_edges.get(&e.key()))
      .map(|e| e.fitch_subs().len())
      .sum();
    assert_eq!(total_subs, 4, "mutation content is unchanged");

    Ok(())
  }
}
