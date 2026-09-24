#[cfg(test)]
mod tests {
  use crate::optimize::topology::hoist_reversions::{hoist_reverting_child, slide_bifurcating_root_for_child};
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::seq::indel::InDel;
  use crate::seq::indel::InDelKind;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  use helpers::{Hoisted, c, edge_indels, edge_subs, make_partition, sub};

  const NWK: &str = "(((A:0.1,B:0.1,Z:0.1)V:0.2)U:0.1)root:0.0;";

  #[test]
  fn test_hoist_reversions_large_t_not_duplicated() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let va = find_edge_key(&graph, &names, "V", "A").unwrap();

    let (partition, _) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        (
          "U",
          "V",
          vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G'), sub(b'G', 10, b'A')],
        ),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let mut sparse = vec![partition];

    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, uv, va, &mut branch_lengths)?;

    let h = Hoisted::locate(&graph, &names, "V", "A");
    let p = &sparse[0];
    assert_eq!(edge_subs(p, h.un), vec![sub(b'C', 5, b'G'), sub(b'G', 10, b'A')]);
    assert_eq!(edge_subs(p, h.nv), vec![sub(b'A', 0, b'T')]);
    assert_eq!(edge_subs(p, h.nc), Vec::<Sub>::new());

    assert_eq!(graph.get_node(h.v).unwrap().degree_out(), 2);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_chain_composed() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let va = find_edge_key(&graph, &names, "V", "A").unwrap();

    let (partition, _) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "A", vec![sub(b'T', 0, b'G')]),
      ],
    );
    let mut sparse = vec![partition];

    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, uv, va, &mut branch_lengths)?;

    let h = Hoisted::locate(&graph, &names, "V", "A");
    let p = &sparse[0];
    assert_eq!(edge_subs(p, h.un), Vec::<Sub>::new());
    assert_eq!(edge_subs(p, h.nv), vec![sub(b'A', 0, b'T')]);
    assert_eq!(edge_subs(p, h.nc), vec![sub(b'A', 0, b'G')]);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_reversion_removed_reduces_count() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let va = find_edge_key(&graph, &names, "V", "A").unwrap();

    let (partition, _) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let mut sparse = vec![partition];

    let before = helpers::total_subs(&graph, &sparse[0]);
    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, uv, va, &mut branch_lengths)?;
    let after = helpers::total_subs(&graph, &sparse[0]);

    assert_eq!(before, 2);
    assert_eq!(after, 1);

    let h = Hoisted::locate(&graph, &names, "V", "A");
    let p = &sparse[0];
    assert_eq!(edge_subs(p, h.nv), vec![sub(b'A', 0, b'T')]);
    assert_eq!(edge_subs(p, h.nc), Vec::<Sub>::new());
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_branch_length_distance_preserved() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let va = find_edge_key(&graph, &names, "V", "A").unwrap();
    let ru = find_edge_key(&graph, &names, "root", "U").unwrap();

    let (partition, _) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        (
          "U",
          "V",
          vec![sub(b'A', 0, b'T'), sub(b'C', 5, b'G'), sub(b'G', 10, b'A')],
        ),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let mut sparse = vec![partition];

    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, uv, va, &mut branch_lengths)?;

    let h = Hoisted::locate(&graph, &names, "V", "A");
    let bl = |ek: GraphEdgeKey| branch_lengths[&ek].unwrap_or(0.0);
    let root_to_v = bl(ru) + bl(h.un) + bl(h.nv);
    let root_to_a = bl(ru) + bl(h.un) + bl(h.nc);

    assert_abs_diff_eq!(root_to_v, 0.1 + 0.2, epsilon = 1e-9);
    assert_abs_diff_eq!(root_to_a, 0.1 + 0.2 + 0.1, epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_multi_partition() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let va = find_edge_key(&graph, &names, "V", "A").unwrap();

    let (p0, _) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T'), sub(b'G', 10, b'C')]),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let (p1, _) = make_partition(
      &graph,
      &names,
      1,
      100,
      &[
        ("U", "V", vec![sub(b'C', 5, b'G')]),
        ("V", "A", vec![sub(b'G', 5, b'C')]),
      ],
    );
    let mut sparse = vec![p0, p1];

    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, uv, va, &mut branch_lengths)?;

    let h = Hoisted::locate(&graph, &names, "V", "A");
    let g0 = &sparse[0];
    assert_eq!(edge_subs(g0, h.un), vec![sub(b'G', 10, b'C')]);
    assert_eq!(edge_subs(g0, h.nv), vec![sub(b'A', 0, b'T')]);
    assert_eq!(edge_subs(g0, h.nc), Vec::<Sub>::new());

    let g1 = &sparse[1];
    assert_eq!(edge_subs(g1, h.un), Vec::<Sub>::new());
    assert_eq!(edge_subs(g1, h.nv), vec![sub(b'C', 5, b'G')]);
    assert_eq!(edge_subs(g1, h.nc), Vec::<Sub>::new());
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_indel_cancellation() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let va = find_edge_key(&graph, &names, "V", "A").unwrap();

    let (mut partition, _) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let del = InDel::del((20, 23), [c(b'A'), c(b'A'), c(b'A')].as_slice())?;
    let ins = InDel::ins((20, 23), [c(b'A'), c(b'A'), c(b'A')].as_slice())?;
    {
      let p = &mut partition;
      p.obs_edges.get_mut(&uv).unwrap().indels = vec![del.clone()];
      p.obs_edges.get_mut(&va).unwrap().indels = vec![ins];
    }
    let mut sparse = vec![partition];

    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, uv, va, &mut branch_lengths)?;

    let h = Hoisted::locate(&graph, &names, "V", "A");
    let p = &sparse[0];
    assert_eq!(edge_indels(p, h.un), Vec::<InDel>::new());
    assert_eq!(edge_indels(p, h.nv), vec![del]);
    assert_eq!(edge_indels(p, h.nc), Vec::<InDel>::new());
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_indel_overlap_fallback() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let va = find_edge_key(&graph, &names, "V", "A").unwrap();

    let (mut partition, _) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let parent_del = InDel::del((20, 25), [c(b'A'); 5].as_slice())?;
    let child_del = InDel::del((22, 28), [c(b'A'); 6].as_slice())?;
    {
      let p = &mut partition;
      p.obs_edges.get_mut(&uv).unwrap().indels = vec![parent_del.clone()];
      p.obs_edges.get_mut(&va).unwrap().indels = vec![child_del];
    }
    let mut sparse = vec![partition];

    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, uv, va, &mut branch_lengths)?;

    let h = Hoisted::locate(&graph, &names, "V", "A");
    let p = &sparse[0];
    assert_eq!(edge_indels(p, h.un), Vec::<InDel>::new());
    assert_eq!(edge_indels(p, h.nv), vec![parent_del]);

    let nc = edge_indels(p, h.nc);
    assert_eq!(nc.len(), 1);
    assert_eq!(nc[0].range, (20, 28));
    assert_eq!(nc[0].kind, InDelKind::Deletion);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_indel_no_interaction_hoisted() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let uv = find_edge_key(&graph, &names, "U", "V").unwrap();
    let va = find_edge_key(&graph, &names, "V", "A").unwrap();

    let (mut partition, _) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let parent_del = InDel::del((20, 23), [c(b'A'); 3].as_slice())?;
    let child_del = InDel::del((50, 53), [c(b'A'); 3].as_slice())?;
    {
      let p = &mut partition;
      p.obs_edges.get_mut(&uv).unwrap().indels = vec![parent_del.clone()];
      p.obs_edges.get_mut(&va).unwrap().indels = vec![child_del.clone()];
    }
    let mut sparse = vec![partition];

    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, uv, va, &mut branch_lengths)?;

    let h = Hoisted::locate(&graph, &names, "V", "A");
    let p = &sparse[0];
    assert_eq!(edge_indels(p, h.un), vec![parent_del]);
    assert_eq!(edge_indels(p, h.nv), Vec::<InDel>::new());
    assert_eq!(edge_indels(p, h.nc), vec![child_del]);
    Ok(())
  }

  const NWK_BIFURCATING: &str = "((C1:0.1,C2:0.1)V:0.1,S:0.1)root:0.0;";

  #[test]
  fn test_hoist_reversions_slide_moves_sibling_sub_to_parent() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK_BIFURCATING)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_v = find_edge_key(&graph, &names, "root", "V").unwrap();
    let root_s = find_edge_key(&graph, &names, "root", "S").unwrap();
    let v_c1 = find_edge_key(&graph, &names, "V", "C1").unwrap();

    let (partition, node_states) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("root", "S", vec![sub(b'A', 3, b'G')]),
        ("V", "C1", vec![sub(b'A', 3, b'G')]),
      ],
    );
    let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

    slide_bifurcating_root_for_child(&mut sparse, &mut node_states, root_key, root_v, root_s, v_c1)?;

    let p = &sparse[0];
    assert_eq!(p.root_sequence[3], c(b'G'));
    assert_eq!(node_states[0][&root_key].sequence[3], c(b'G'));
    assert_eq!(edge_subs(p, root_s), Vec::<Sub>::new());
    assert_eq!(edge_subs(p, root_v), vec![sub(b'G', 3, b'A')]);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_slide_then_hoist_removes_cross_root_reversion() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(NWK_BIFURCATING)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let root_key = find_node_key_by_name(&graph, &names, "root").unwrap();
    let root_v = find_edge_key(&graph, &names, "root", "V").unwrap();
    let root_s = find_edge_key(&graph, &names, "root", "S").unwrap();
    let v_c1 = find_edge_key(&graph, &names, "V", "C1").unwrap();

    let (partition, node_states) = make_partition(
      &graph,
      &names,
      0,
      100,
      &[
        ("root", "S", vec![sub(b'A', 3, b'G')]),
        ("V", "C1", vec![sub(b'A', 3, b'G')]),
      ],
    );
    let mut sparse = vec![partition];
    let mut node_states = vec![node_states];

    let before = helpers::total_subs(&graph, &sparse[0]);
    slide_bifurcating_root_for_child(&mut sparse, &mut node_states, root_key, root_v, root_s, v_c1)?;
    let after_slide = helpers::total_subs(&graph, &sparse[0]);
    let mut branch_lengths = branch_lengths;
    hoist_reverting_child(&mut graph, &mut sparse, root_v, v_c1, &mut branch_lengths)?;
    let after_hoist = helpers::total_subs(&graph, &sparse[0]);

    assert_eq!(before, 2);
    assert_eq!(after_slide, 2);
    assert_eq!(after_hoist, 1);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};

    use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
    use maplit::btreemap;
    use treetime_primitives::{AsciiChar, Seq};

    pub(super) fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub(super) fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(c(reff), pos, c(qry)).unwrap()
    }

    pub(super) struct Hoisted {
      pub v: GraphNodeKey,
      pub un: GraphEdgeKey,
      pub nv: GraphEdgeKey,
      pub nc: GraphEdgeKey,
    }

    impl Hoisted {
      pub(crate) fn locate(
        graph: &Graph,
        names: &BTreeMap<GraphNodeKey, Option<String>>,
        v_name: &str,
        c_name: &str,
      ) -> Self {
        let v = find_node_key_by_name(graph, names, v_name).unwrap();
        let c = find_node_key_by_name(graph, names, c_name).unwrap();
        let nv = single_inbound(graph, v);
        let n = graph.get_source_node_key(nv).unwrap();
        let un = single_inbound(graph, n);
        let nc = single_inbound(graph, c);
        Self { v, un, nv, nc }
      }
    }

    fn single_inbound(graph: &Graph, node_key: GraphNodeKey) -> GraphEdgeKey {
      let node = graph.get_node(node_key).unwrap();
      match node.inbound() {
        [edge_key] => *edge_key,
        other => panic!("expected exactly one inbound edge, found {}", other.len()),
      }
    }

    pub(super) fn edge_subs(recon: &PartitionMarginalSparse, edge_key: GraphEdgeKey) -> Vec<Sub> {
      recon.obs_edges[&edge_key].fitch_subs().to_vec()
    }

    pub(super) fn edge_indels(recon: &PartitionMarginalSparse, edge_key: GraphEdgeKey) -> Vec<InDel> {
      recon.obs_edges[&edge_key].indels.clone()
    }

    pub(super) fn total_subs(graph: &Graph, recon: &PartitionMarginalSparse) -> usize {
      graph
        .get_edges()
        .filter_map(|e| recon.obs_edges.get(&e.key()))
        .map(|e| e.fitch_subs().len())
        .sum()
    }

    pub(super) fn make_partition(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      index: usize,
      length: usize,
      edge_mutations: &[(&str, &str, Vec<Sub>)],
    ) -> (PartitionMarginalSparse, BTreeMap<GraphNodeKey, SparseNodeState>) {
      let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();

      let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
      for (_, _, subs) in edge_mutations {
        for s in subs {
          if s.pos() < length {
            ref_seq[s.pos()] = s.reff();
          }
        }
      }

      let mut obs_nodes = btreemap! {};
      let mut node_states = btreemap! {};
      for node in graph.get_nodes() {
        let key = node.key();
        obs_nodes.insert(key, SparseNodeObs::empty(&alphabet));
        node_states.insert(key, SparseNodeState::leaf(&ref_seq));
      }

      let mut obs_edges = btreemap! {};
      for (source, target, subs) in edge_mutations {
        let edge_key =
          find_edge_key(graph, names, source, target).unwrap_or_else(|| panic!("edge {source}->{target} missing"));
        obs_edges.insert(edge_key, SparseEdgeObs::with_fitch_subs(subs.clone()));
      }

      let partition = PartitionMarginalSparse {
        index,
        alphabet,
        length,
        root_sequence: ref_seq,
        obs_nodes,
        obs_edges,
      };

      (partition, node_states)
    }
  }
}
