#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};

  use crate::optimize::topology::collapse::collapse_edge;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs};
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AsciiChar;
  use treetime_primitives::seq;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  fn populate_test_nodes(recon: &mut PartitionMarginalSparse, graph: &Graph) {
    let ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A')).take(recon.length).collect();
    if recon.root_sequence.is_empty() {
      recon.root_sequence = ref_seq;
    }
    let alphabet = recon.alphabet.clone();
    for node in graph.get_nodes() {
      let key = node.key();
      recon
        .obs_nodes
        .entry(key)
        .or_insert_with(|| SparseNodeObs::empty(&alphabet));
    }
  }

  fn make_sparse_reconstruction(length: usize) -> Result<PartitionMarginalSparse, Report> {
    let partition = PartitionMarginalSparse {
      index: 0,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length,
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
      root_sequence: seq![],
    };
    Ok(partition)
  }

  #[test]
  fn test_topology_collapse_edge_sparse_composes_subs() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();

    let mut recon = make_sparse_reconstruction(100)?;
    populate_test_nodes(&mut recon, &graph);

    recon
      .obs_edges
      .insert(ri_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    recon
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));
    recon.obs_edges.insert(ib_key, SparseEdgeObs::default());

    let mut sparse = vec![recon];

    let i_node_key = find_node_key_by_name(&graph, &names, "I").unwrap();

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    let p = &sparse[0];
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.key()).cloned().flatten();
      let edge_data = &p.obs_edges[&edge.key()];
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

    assert!(!p.obs_nodes.contains_key(&i_node_key));
    assert!(!p.obs_edges.contains_key(&ri_key));

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_graph_cleanup() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let i_key = find_node_key_by_name(&graph, &names, "I").unwrap();

    let mut sparse: Vec<PartitionMarginalSparse> = vec![];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    assert!(
      graph.get_node(i_key).is_none(),
      "removed node should be gone from graph"
    );
    assert!(
      graph.get_edge(ri_key).is_none(),
      "removed edge should be gone from graph"
    );

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_branch_length_sum() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut sparse: Vec<PartitionMarginalSparse> = vec![];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.key()).cloned().flatten();
      let bl = branch_lengths[&edge.key()];

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
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut sparse: Vec<PartitionMarginalSparse> = vec![];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.key()).cloned().flatten();
      let bl = branch_lengths[&edge.key()];

      match target_name.as_deref() {
        Some("A") => assert_abs_diff_eq!(bl.unwrap(), 0.1, epsilon = 1e-7),
        Some("B") => assert_abs_diff_eq!(bl.unwrap(), 0.2, epsilon = 1e-7),
        other => unreachable!("unexpected target node: {other:?}"),
      }
    }

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_branch_length_none_plus_some() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();

    let mut sparse: Vec<PartitionMarginalSparse> = vec![];

    let mut branch_lengths = branch_lengths;
    branch_lengths.insert(ia_key, None);
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.key()).cloned().flatten();
      match target_name.as_deref() {
        Some("A") => assert_eq!(branch_lengths[&edge.key()], None),
        Some("B") => assert_abs_diff_eq!(branch_lengths[&edge.key()].unwrap(), 0.5, epsilon = 1e-7),
        other => unreachable!("unexpected target node: {other:?}"),
      }
    }

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_branch_length_some_plus_none() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut sparse: Vec<PartitionMarginalSparse> = vec![];

    let mut branch_lengths = branch_lengths;
    branch_lengths.insert(ri_key, None);
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.key()).cloned().flatten();
      match target_name.as_deref() {
        Some("A") => assert_abs_diff_eq!(branch_lengths[&edge.key()].unwrap(), 0.1, epsilon = 1e-7),
        Some("B") => assert_abs_diff_eq!(branch_lengths[&edge.key()].unwrap(), 0.2, epsilon = 1e-7),
        other => unreachable!("unexpected target node: {other:?}"),
      }
    }

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_indel_concatenation() -> Result<(), Report> {
    use crate::seq::indel::InDel;

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();

    let mut recon = make_sparse_reconstruction(100)?;
    populate_test_nodes(&mut recon, &graph);

    let collapsed_indel = InDel::ins((0, 3), [c(b'A'), c(b'C'), c(b'G')].as_slice()).unwrap();
    let child_a_indel = InDel::del((10, 12), [c(b'T'), c(b'T')].as_slice()).unwrap();

    recon.obs_edges.insert(
      ri_key,
      SparseEdgeObs::with_fitch_subs_and_indels(vec![], vec![collapsed_indel.clone()]),
    );
    recon.obs_edges.insert(
      ia_key,
      SparseEdgeObs::with_fitch_subs_and_indels(vec![], vec![child_a_indel.clone()]),
    );
    recon.obs_edges.insert(ib_key, SparseEdgeObs::default());

    let mut sparse = vec![recon];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    let p = &sparse[0];
    for edge in graph.get_edges() {
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.key()).cloned().flatten();
      let edge_data = &p.obs_edges[&edge.key()];
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
    let nwk_parsed = nwk_read_str("((A:0.1)I:0.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();

    let mut recon = make_sparse_reconstruction(100)?;
    populate_test_nodes(&mut recon, &graph);

    recon
      .obs_edges
      .insert(ri_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    recon
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));

    let mut sparse = vec![recon];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    let p = &sparse[0];
    let root_to_a_edge = graph
      .get_edges()
      .find(|e| {
        let t = e.target();
        graph
          .get_node(t)
          .and_then(|n| names.get(&n.key()).cloned().flatten())
          .as_deref()
          == Some("A")
      })
      .cloned()
      .unwrap();
    let edge_data = &p.obs_edges[&root_to_a_edge.key()];
    let expected: &[Sub] = &[];
    assert_eq!(edge_data.fitch_subs(), expected);

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_no_partitions() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let i_node_key = find_node_key_by_name(&graph, &names, "I").unwrap();

    let mut sparse: Vec<PartitionMarginalSparse> = vec![];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    assert!(graph.get_node(i_node_key).is_none());
    assert_eq!(graph.get_nodes().count(), 3);
    assert_eq!(graph.get_edges().count(), 2);

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_multiple_sparse_partitions() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();

    let mut recon_a = make_sparse_reconstruction(100)?;
    populate_test_nodes(&mut recon_a, &graph);
    recon_a
      .obs_edges
      .insert(ri_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    recon_a.obs_edges.insert(ia_key, SparseEdgeObs::default());
    recon_a.obs_edges.insert(ib_key, SparseEdgeObs::default());

    let mut recon_b = make_sparse_reconstruction(100)?;
    recon_b.index = 1;
    populate_test_nodes(&mut recon_b, &graph);
    recon_b
      .obs_edges
      .insert(ri_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));
    recon_b.obs_edges.insert(ia_key, SparseEdgeObs::default());
    recon_b.obs_edges.insert(ib_key, SparseEdgeObs::default());

    let mut sparse = vec![recon_a, recon_b];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    let p0 = &sparse[0];
    for edge in graph.get_edges() {
      let data = &p0.obs_edges[&edge.key()];
      assert_eq!(
        data.fitch_subs(),
        &[sub(b'A', 0, b'T')],
        "partition 0: each child inherits A0T"
      );
    }

    let p1 = &sparse[1];
    for edge in graph.get_edges() {
      let data = &p1.obs_edges[&edge.key()];
      assert_eq!(
        data.fitch_subs(),
        &[sub(b'G', 5, b'C')],
        "partition 1: each child inherits G5C"
      );
    }

    Ok(())
  }
}
