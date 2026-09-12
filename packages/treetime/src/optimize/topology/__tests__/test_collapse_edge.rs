#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::topology::collapse::collapse_edge;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs};
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_primitives::AsciiChar;
  use treetime_primitives::seq;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
    Sub::new(c(reff), pos, c(qry)).unwrap()
  }

  fn populate_test_nodes(recon: &mut SparseReconstruction, graph: &Graph) {
    let ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A')).take(recon.partition.length).collect();
    if recon.partition.root_sequence.is_empty() {
      recon.partition.root_sequence = ref_seq.clone();
    }
    let alphabet = recon.partition.alphabet.clone();
    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      recon
        .partition
        .obs_nodes
        .entry(key)
        .or_insert_with(|| SparseNodeObs::empty(&alphabet));
    }
  }

  fn make_sparse_reconstruction(length: usize) -> Result<SparseReconstruction, Report> {
    let partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: Alphabet::new(AlphabetName::Nuc)?,
      length,
      obs_nodes: btreemap! {},
      obs_edges: btreemap! {},
      root_sequence: seq![],
    };
    Ok(SparseReconstruction {
      partition,
      node_states: BTreeMap::new(),
      backward: BTreeMap::new(),
      forward: BTreeMap::new(),
      estimates: BTreeMap::new(),
    })
  }

  #[test]
  fn test_topology_collapse_edge_sparse_composes_subs() -> Result<(), Report> {
    // Tree: root -> I (bl=0.0) -> A, B
    // I has sub A0T; A has sub G5C; B has no subs
    // After collapse: root -> A has {A0T, G5C}, root -> B has {A0T}
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();

    let mut recon = make_sparse_reconstruction(100)?;
    populate_test_nodes(&mut recon, &graph);

    recon
      .partition
      .obs_edges
      .insert(ri_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    recon
      .partition
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));
    recon.partition.obs_edges.insert(ib_key, SparseEdgeObs::default());

    let mut sparse = vec![recon];

    let i_node_key = find_node_key_by_name(&graph, &names, "I").unwrap();

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2);

    let p = &sparse[0].partition;
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.read_arc().key()).cloned().flatten();
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

    // Stale observations removed
    assert!(!p.obs_nodes.contains_key(&i_node_key));
    assert!(!p.obs_edges.contains_key(&ri_key));

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_graph_cleanup() -> Result<(), Report> {
    // Graph-level cleanup: the target node and collapsed edge are removed from the graph. Dense
    // per-node/per-edge state is reconciled centrally after a topology batch (reconcile_dense_family
    // in run_loop), not by collapse_edge, so this checks only the graph effect collapse_edge owns.
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let i_key = find_node_key_by_name(&graph, &names, "I").unwrap();

    let mut sparse: Vec<SparseReconstruction> = vec![];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    assert!(graph.get_node(i_key).is_none(), "removed node should be gone from graph");
    assert!(graph.get_edge(ri_key).is_none(), "removed edge should be gone from graph");

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_branch_length_sum() -> Result<(), Report> {
    // Collapsed edge has bl=0.3, child edges bl=0.1 and bl=0.2
    // After collapse: child edges bl = 0.3 + 0.1 = 0.4 and 0.3 + 0.2 = 0.5
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut sparse: Vec<SparseReconstruction> = vec![];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.read_arc().key()).cloned().flatten();
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
    // Collapsed edge has bl=0.0, child edges preserved unchanged.
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)I:0.0)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut sparse: Vec<SparseReconstruction> = vec![];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.read_arc().key()).cloned().flatten();
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
    // Collapsed edge length present, one child length missing (None): the missing child length
    // is preserved as None (no sum), the present child is summed. Oracle: the Option-aware sum
    // in `collapse_edge` only sums when both operands are `Some`.
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let mut graph: Graph = graph;
    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();

    let mut sparse: Vec<SparseReconstruction> = vec![];

    let mut branch_lengths = branch_lengths;
    branch_lengths.insert(ia_key, None);
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.read_arc().key()).cloned().flatten();
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
    // Collapsed edge length missing (None): child lengths are preserved unchanged. Oracle: the
    // Option-aware sum in `collapse_edge` only sums when both operands are `Some`.
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let mut graph: Graph = graph;
    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();

    let mut sparse: Vec<SparseReconstruction> = vec![];

    let mut branch_lengths = branch_lengths;
    branch_lengths.insert(ri_key, None);
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.read_arc().key()).cloned().flatten();
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
    // Indels on the collapsed edge must be preserved on each former-child edge
    // (collapsed-edge indels prepended to child indels).
    use crate::seq::indel::InDel;

    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;

    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();

    let mut recon = make_sparse_reconstruction(100)?;
    populate_test_nodes(&mut recon, &graph);

    let collapsed_indel = InDel::ins((0, 3), [c(b'A'), c(b'C'), c(b'G')].as_slice()).unwrap();
    let child_a_indel = InDel::del((10, 12), [c(b'T'), c(b'T')].as_slice()).unwrap();

    recon.partition.obs_edges.insert(
      ri_key,
      SparseEdgeObs::with_fitch_subs_and_indels(vec![], vec![collapsed_indel.clone()]),
    );
    recon.partition.obs_edges.insert(
      ia_key,
      SparseEdgeObs::with_fitch_subs_and_indels(vec![], vec![child_a_indel.clone()]),
    );
    recon.partition.obs_edges.insert(ib_key, SparseEdgeObs::default());

    let mut sparse = vec![recon];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    let p = &sparse[0].partition;
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let target = graph.get_node(edge.target()).unwrap();
      let target_name = names.get(&target.read_arc().key()).cloned().flatten();
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
    // Collapsed edge A0T + child edge T0A = no net change (reversion).
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1)I:0.0)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();

    let mut recon = make_sparse_reconstruction(100)?;
    populate_test_nodes(&mut recon, &graph);

    recon
      .partition
      .obs_edges
      .insert(ri_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    recon
      .partition
      .obs_edges
      .insert(ia_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'T', 0, b'A')]));

    let mut sparse = vec![recon];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    let p = &sparse[0].partition;
    let root_to_a_edge = graph
      .get_edges()
      .iter()
      .find(|e| {
        let t = e.read_arc().target();
        graph
          .get_node(t)
          .and_then(|n| names.get(&n.read_arc().key()).cloned().flatten())
          .as_deref()
          == Some("A")
      })
      .cloned()
      .unwrap();
    let edge_data = &p.obs_edges[&root_to_a_edge.read_arc().key()];
    let expected: &[Sub] = &[];
    assert_eq!(edge_data.fitch_subs(), expected);

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_no_partitions() -> Result<(), Report> {
    // Graph-only collapse with no partitions: topology still changes correctly.
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.2)I:0.3)root;")?;
    let mut graph: Graph = graph;
    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let i_node_key = find_node_key_by_name(&graph, &names, "I").unwrap();

    let mut sparse: Vec<SparseReconstruction> = vec![];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    assert!(graph.get_node(i_node_key).is_none());
    assert_eq!(graph.get_nodes().len(), 3); // root, A, B
    assert_eq!(graph.get_edges().len(), 2);

    Ok(())
  }

  #[test]
  fn test_topology_collapse_edge_multiple_sparse_partitions() -> Result<(), Report> {
    // Two sparse partitions with independent edge data should both be updated.
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.1)I:0.0)root;")?;
    let mut graph: Graph = graph;

    let ri_key = find_edge_key(&graph, &names, "root", "I").unwrap();
    let ia_key = find_edge_key(&graph, &names, "I", "A").unwrap();
    let ib_key = find_edge_key(&graph, &names, "I", "B").unwrap();

    let mut recon_a = make_sparse_reconstruction(100)?;
    populate_test_nodes(&mut recon_a, &graph);
    recon_a
      .partition
      .obs_edges
      .insert(ri_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'A', 0, b'T')]));
    recon_a.partition.obs_edges.insert(ia_key, SparseEdgeObs::default());
    recon_a.partition.obs_edges.insert(ib_key, SparseEdgeObs::default());

    let mut recon_b = make_sparse_reconstruction(100)?;
    recon_b.partition.index = 1;
    populate_test_nodes(&mut recon_b, &graph);
    recon_b
      .partition
      .obs_edges
      .insert(ri_key, SparseEdgeObs::with_fitch_subs(vec![sub(b'G', 5, b'C')]));
    recon_b.partition.obs_edges.insert(ia_key, SparseEdgeObs::default());
    recon_b.partition.obs_edges.insert(ib_key, SparseEdgeObs::default());

    let mut sparse = vec![recon_a, recon_b];

    let mut branch_lengths = branch_lengths;
    collapse_edge(&mut graph, &mut sparse, ri_key, &mut branch_lengths)?;
    graph.build()?;

    let p0 = &sparse[0].partition;
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let data = &p0.obs_edges[&edge.key()];
      assert_eq!(
        data.fitch_subs(),
        &[sub(b'A', 0, b'T')],
        "partition 0: each child inherits A0T"
      );
    }

    let p1 = &sparse[1].partition;
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
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
