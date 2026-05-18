#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::topology::collapse::collapse_edge;
  use crate::partition::marginal_dense::PartitionMarginalDense;

  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::payload::ancestral::GraphAncestral;
  use crate::partition::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
  use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AsciiChar;
  use treetime_primitives::seq;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  fn populate_test_nodes(partition: &mut PartitionMarginalSparse, graph: &GraphAncestral) {
    let ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A')).take(partition.length).collect();
    if partition.root_sequence.is_empty() {
      partition.root_sequence = ref_seq.clone();
    }
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      partition.nodes.entry(key).or_insert_with(|| {
        let mut node_part = SparseNodePartition::empty(&partition.alphabet);
        node_part.seq.sequence = ref_seq.clone();
        node_part
      });
    }
  }

  fn make_sparse_partition(length: usize) -> Result<PartitionMarginalSparse, Report> {
    Ok(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    })
  }

  #[test]
  fn test_topology_collapse_edge_sparse_composes_subs() -> Result<(), Report> {
    // Tree: root -> I (bl=0.0) -> A, B
    // I has sub A0T; A has sub G5C; B has no subs
    // After collapse: root -> A has {A0T, G5C}, root -> B has {A0T}
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, "I", "B").unwrap();

    let mut partition = make_sparse_partition(100)?;
    populate_test_nodes(&mut partition, &graph);

    partition
      .edges
      .insert(ri_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .edges
      .insert(ia_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'G', 5, b'C')]));
    partition.edges.insert(ib_key, SparseEdgePartition::default());

    let sparse = vec![Arc::new(RwLock::new(partition))];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    let i_node_key = find_node_key_by_name(&graph, "I").unwrap();

    collapse_edge(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2);

    let p = sparse[0].read_arc();
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let edge_data = &p.edges[&edge.key()];
      match target_name.as_deref() {
        Some("A") => {
          assert_eq!(edge_data.fitch_subs(), &[sub(b'A', 0, b'T'), sub(b'G', 5, b'C')]);
        },
        Some("B") => {
          assert_eq!(edge_data.fitch_subs(), &[sub(b'A', 0, b'T')]);
        },
        other => unreachable!("unexpected target node: {other:?}"),
      }
    }

    // Stale entries removed
    assert!(!p.nodes.contains_key(&i_node_key));
    assert!(!p.edges.contains_key(&ri_key));

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_dense_cleanup() -> Result<(), Report> {
    // Dense partition: stale node/edge entries should be removed after collapse
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();
    let i_key = find_node_key_by_name(&graph, "I").unwrap();

    let mut dense_partition = PartitionMarginalDense::new(0, jc69(JC69Params::default())?, Alphabet::new(AlphabetName::Nuc)?, 10);

    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      dense_partition.data.nodes.insert(
        key,
        DenseNodePartition {
          seq: DenseSeqInfo::default(),
          profile: DenseSeqDistribution::default(),
        },
      );
    }
    for edge in graph.get_edges() {
      let key = edge.read_arc().key();
      dense_partition.data.edges.insert(key, DenseEdgePartition::default());
    }

    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let dense = vec![Arc::new(RwLock::new(dense_partition))];

    collapse_edge(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    let d = dense[0].read_arc();
    assert!(!d.data.nodes.contains_key(&i_key), "removed node should be cleaned up");
    assert!(!d.data.edges.contains_key(&ri_key), "removed edge should be cleaned up");

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_branch_length_sum() -> Result<(), Report> {
    // Collapsed edge has bl=0.3, child edges bl=0.1 and bl=0.2
    // After collapse: child edges bl = 0.3 + 0.1 = 0.4 and 0.3 + 0.2 = 0.5
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();

    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    collapse_edge(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length();

      match target_name.as_deref() {
        Some("A") => assert_abs_diff_eq!(bl.unwrap(), 0.4, epsilon = 1e-7),
        Some("B") => assert_abs_diff_eq!(bl.unwrap(), 0.5, epsilon = 1e-7),
        other => unreachable!("unexpected target node: {other:?}"),
      }
    }

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_branch_length_sum_with_zero() -> Result<(), Report> {
    // Collapsed edge has bl=0.0, child edges preserved unchanged.
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();

    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    collapse_edge(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let bl = edge.payload().read_arc().branch_length();

      match target_name.as_deref() {
        Some("A") => assert_abs_diff_eq!(bl.unwrap(), 0.1, epsilon = 1e-7),
        Some("B") => assert_abs_diff_eq!(bl.unwrap(), 0.2, epsilon = 1e-7),
        other => unreachable!("unexpected target node: {other:?}"),
      }
    }

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_indel_concatenation() -> Result<(), Report> {
    // Indels on the collapsed edge must be preserved on each former-child edge
    // (collapsed-edge indels prepended to child indels).
    use crate::seq::indel::InDel;

    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, "I", "B").unwrap();

    let mut partition = make_sparse_partition(100)?;
    populate_test_nodes(&mut partition, &graph);

    let collapsed_indel = InDel::ins((0, 3), [c(b'A'), c(b'C'), c(b'G')].as_slice());
    let child_a_indel = InDel::del((10, 12), [c(b'T'), c(b'T')].as_slice());

    partition.edges.insert(
      ri_key,
      SparseEdgePartition::with_fitch_subs_and_indels(vec![], vec![collapsed_indel.clone()]),
    );
    partition.edges.insert(
      ia_key,
      SparseEdgePartition::with_fitch_subs_and_indels(vec![], vec![child_a_indel.clone()]),
    );
    partition.edges.insert(ib_key, SparseEdgePartition::default());

    let sparse = vec![Arc::new(RwLock::new(partition))];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    collapse_edge(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    let p = sparse[0].read_arc();
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = target.read_arc().payload().read_arc().name.clone();
      let edge_data = &p.edges[&edge.key()];
      match target_name.as_deref() {
        Some("A") => {
          assert_eq!(edge_data.indels, vec![collapsed_indel.clone(), child_a_indel.clone()]);
        },
        Some("B") => {
          assert_eq!(edge_data.indels, vec![collapsed_indel.clone()]);
        },
        other => unreachable!("unexpected target node: {other:?}"),
      }
    }

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_reversion_cancels() -> Result<(), Report> {
    // Collapsed edge A0T + child edge T0A = no net change (reversion).
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, "I", "A").unwrap();

    let mut partition = make_sparse_partition(100)?;
    populate_test_nodes(&mut partition, &graph);

    partition
      .edges
      .insert(ri_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition
      .edges
      .insert(ia_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'T', 0, b'A')]));

    let sparse = vec![Arc::new(RwLock::new(partition))];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    collapse_edge(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    let p = sparse[0].read_arc();
    let root_to_a_edge = graph
      .get_edges()
      .iter()
      .find(|e| {
        let t = e.read_arc().target();
        graph
          .get_node(t)
          .and_then(|n| n.read_arc().payload().read_arc().name.clone())
          .as_deref()
          == Some("A")
      })
      .cloned()
      .unwrap();
    let edge_data = &p.edges[&root_to_a_edge.read_arc().key()];
    let expected: &[Sub] = &[];
    assert_eq!(edge_data.fitch_subs(), expected);

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_no_partitions() -> Result<(), Report> {
    // Graph-only collapse with no partitions: topology still changes correctly.
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let ri_key = find_edge_key(&graph, "root", "I").unwrap();
    let i_node_key = find_node_key_by_name(&graph, "I").unwrap();

    let sparse: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    collapse_edge(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    assert!(graph.get_node(i_node_key).is_none());
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2);

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_multiple_sparse_partitions() -> Result<(), Report> {
    // Two sparse partitions with independent edge data should both be updated.
    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;

    let ri_key = find_edge_key(&graph, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, "I", "B").unwrap();

    let mut partition_a = make_sparse_partition(100)?;
    populate_test_nodes(&mut partition_a, &graph);
    partition_a
      .edges
      .insert(ri_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    partition_a.edges.insert(ia_key, SparseEdgePartition::default());
    partition_a.edges.insert(ib_key, SparseEdgePartition::default());

    let mut partition_b = make_sparse_partition(100)?;
    partition_b.index = 1;
    populate_test_nodes(&mut partition_b, &graph);
    partition_b
      .edges
      .insert(ri_key, SparseEdgePartition::with_fitch_subs(vec![sub(b'G', 5, b'C')]));
    partition_b.edges.insert(ia_key, SparseEdgePartition::default());
    partition_b.edges.insert(ib_key, SparseEdgePartition::default());

    let sparse = vec![Arc::new(RwLock::new(partition_a)), Arc::new(RwLock::new(partition_b))];
    let dense: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![];

    collapse_edge(&mut graph, &sparse, &dense, ri_key)?;
    graph.build()?;

    let p0 = sparse[0].read_arc();
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let data = &p0.edges[&edge.key()];
      assert_eq!(
        data.fitch_subs(),
        &[sub(b'A', 0, b'T')],
        "partition 0: each child inherits A0T"
      );
    }
    drop(p0);

    let p1 = sparse[1].read_arc();
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let data = &p1.edges[&edge.key()];
      assert_eq!(
        data.fitch_subs(),
        &[sub(b'G', 5, b'C')],
        "partition 1: each child inherits G5C"
      );
    }

    Ok(())
  }
}
