#[cfg(test)]
mod tests {
  use crate::optimize::topology::hoist_reversions::{hoist_reverting_child, slide_bifurcating_root_for_child};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use crate::seq::indel::InDelKind;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  use helpers::{Hoisted, c, edge_indels, edge_subs, make_partition, no_dense, sub};

  // Tree: root -> U -> V -> {A, B, Z}. The hoist inserts N between U and V, grouping V with A.
  const NWK: &str = "(((A:0.1,B:0.1,Z:0.1)V:0.2)U:0.1)root:0.0;";

  #[test]
  fn test_hoist_reversions_large_t_not_duplicated() -> Result<(), Report> {
    // M_v has three substitutions; the child reverts only one. The two untouched
    // substitutions (T) must land on u->N once and NOT be duplicated onto N->c
    // (which distinguishes the move from re-attaching the child to the parent).
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let uv = find_edge_key(&graph, "U", "V").unwrap();
    let va = find_edge_key(&graph, "V", "A").unwrap();

    let partition = make_partition(
      &graph,
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
    let sparse = vec![partition];

    hoist_reverting_child(&mut graph, &sparse, &no_dense(), uv, va)?;

    let h = Hoisted::locate(&graph, "V", "A");
    let p = sparse[0].read_arc();
    assert_eq!(edge_subs(&p, h.un), vec![sub(b'C', 5, b'G'), sub(b'G', 10, b'A')]);
    assert_eq!(edge_subs(&p, h.nv), vec![sub(b'A', 0, b'T')]);
    assert_eq!(edge_subs(&p, h.nc), Vec::<Sub>::new());

    // V keeps its other children B and Z.
    assert_eq!(graph.get_node(h.v).unwrap().read_arc().degree_out(), 2);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_chain_composed() -> Result<(), Report> {
    // Chain: parent A0T at pos 0, child T0G at pos 0 -> net A0G. The original A0T stays
    // on N->v; the composed A0G moves to N->c.
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let uv = find_edge_key(&graph, "U", "V").unwrap();
    let va = find_edge_key(&graph, "V", "A").unwrap();

    let partition = make_partition(
      &graph,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "A", vec![sub(b'T', 0, b'G')]),
      ],
    );
    let sparse = vec![partition];

    hoist_reverting_child(&mut graph, &sparse, &no_dense(), uv, va)?;

    let h = Hoisted::locate(&graph, "V", "A");
    let p = sparse[0].read_arc();
    assert_eq!(edge_subs(&p, h.un), Vec::<Sub>::new());
    assert_eq!(edge_subs(&p, h.nv), vec![sub(b'A', 0, b'T')]);
    assert_eq!(edge_subs(&p, h.nc), vec![sub(b'A', 0, b'G')]);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_reversion_removed_reduces_count() -> Result<(), Report> {
    // Pure reversion: A0T then T0A. Two mutations before, one after (delta = -1).
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let uv = find_edge_key(&graph, "U", "V").unwrap();
    let va = find_edge_key(&graph, "V", "A").unwrap();

    let partition = make_partition(
      &graph,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T')]),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let sparse = vec![partition];

    let before = helpers::total_subs(&graph, &sparse[0].read_arc());
    hoist_reverting_child(&mut graph, &sparse, &no_dense(), uv, va)?;
    let after = helpers::total_subs(&graph, &sparse[0].read_arc());

    assert_eq!(before, 2);
    assert_eq!(after, 1);

    let h = Hoisted::locate(&graph, "V", "A");
    let p = sparse[0].read_arc();
    assert_eq!(edge_subs(&p, h.nv), vec![sub(b'A', 0, b'T')]);
    assert_eq!(edge_subs(&p, h.nc), Vec::<Sub>::new());
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_branch_length_distance_preserved() -> Result<(), Report> {
    // Distances root->V and root->A are unchanged by the move; the parent edge is split
    // proportionally to substitution count (|T|/|M_v| = 2/3).
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let uv = find_edge_key(&graph, "U", "V").unwrap();
    let va = find_edge_key(&graph, "V", "A").unwrap();
    let ru = find_edge_key(&graph, "root", "U").unwrap();

    let partition = make_partition(
      &graph,
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
    let sparse = vec![partition];

    hoist_reverting_child(&mut graph, &sparse, &no_dense(), uv, va)?;

    let h = Hoisted::locate(&graph, "V", "A");
    let bl = |ek: GraphEdgeKey| helpers::branch_length(&graph, ek);
    let root_to_v = bl(ru) + bl(h.un) + bl(h.nv);
    let root_to_a = bl(ru) + bl(h.un) + bl(h.nc);

    assert_abs_diff_eq!(root_to_v, 0.1 + 0.2, epsilon = 1e-9);
    assert_abs_diff_eq!(root_to_a, 0.1 + 0.2 + 0.1, epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_multi_partition() -> Result<(), Report> {
    // Two partitions revert independent positions. Each partition's edges are split on
    // its own positions; T is per-partition (present in p0, empty in p1).
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let uv = find_edge_key(&graph, "U", "V").unwrap();
    let va = find_edge_key(&graph, "V", "A").unwrap();

    let p0 = make_partition(
      &graph,
      0,
      100,
      &[
        ("U", "V", vec![sub(b'A', 0, b'T'), sub(b'G', 10, b'C')]),
        ("V", "A", vec![sub(b'T', 0, b'A')]),
      ],
    );
    let p1 = make_partition(
      &graph,
      1,
      100,
      &[
        ("U", "V", vec![sub(b'C', 5, b'G')]),
        ("V", "A", vec![sub(b'G', 5, b'C')]),
      ],
    );
    let sparse = vec![p0, p1];

    hoist_reverting_child(&mut graph, &sparse, &no_dense(), uv, va)?;

    let h = Hoisted::locate(&graph, "V", "A");
    let g0 = sparse[0].read_arc();
    assert_eq!(edge_subs(&g0, h.un), vec![sub(b'G', 10, b'C')]);
    assert_eq!(edge_subs(&g0, h.nv), vec![sub(b'A', 0, b'T')]);
    assert_eq!(edge_subs(&g0, h.nc), Vec::<Sub>::new());

    let g1 = sparse[1].read_arc();
    assert_eq!(edge_subs(&g1, h.un), Vec::<Sub>::new());
    assert_eq!(edge_subs(&g1, h.nv), vec![sub(b'C', 5, b'G')]);
    assert_eq!(edge_subs(&g1, h.nc), Vec::<Sub>::new());
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_indel_cancellation() -> Result<(), Report> {
    // A deletion on the parent edge and its inverse insertion on the child edge interact,
    // so the parent indel stays on N->v and the composition cancels on N->c.
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let uv = find_edge_key(&graph, "U", "V").unwrap();
    let va = find_edge_key(&graph, "V", "A").unwrap();

    let partition = make_partition(
      &graph,
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
      let mut p = partition.write_arc();
      p.edges.get_mut(&uv).unwrap().indels = vec![del.clone()];
      p.edges.get_mut(&va).unwrap().indels = vec![ins];
    }
    let sparse = vec![partition];

    hoist_reverting_child(&mut graph, &sparse, &no_dense(), uv, va)?;

    let h = Hoisted::locate(&graph, "V", "A");
    let p = sparse[0].read_arc();
    assert_eq!(edge_indels(&p, h.un), Vec::<InDel>::new());
    assert_eq!(edge_indels(&p, h.nv), vec![del]);
    assert_eq!(edge_indels(&p, h.nc), Vec::<InDel>::new());
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_indel_overlap_fallback() -> Result<(), Report> {
    // Overlapping deletions cannot be cleanly hoisted: the parent deletion stays on N->v
    // and the merged deletion lands on N->c.
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let uv = find_edge_key(&graph, "U", "V").unwrap();
    let va = find_edge_key(&graph, "V", "A").unwrap();

    let partition = make_partition(
      &graph,
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
      let mut p = partition.write_arc();
      p.edges.get_mut(&uv).unwrap().indels = vec![parent_del.clone()];
      p.edges.get_mut(&va).unwrap().indels = vec![child_del];
    }
    let sparse = vec![partition];

    hoist_reverting_child(&mut graph, &sparse, &no_dense(), uv, va)?;

    let h = Hoisted::locate(&graph, "V", "A");
    let p = sparse[0].read_arc();
    assert_eq!(edge_indels(&p, h.un), Vec::<InDel>::new());
    assert_eq!(edge_indels(&p, h.nv), vec![parent_del]);

    let nc = edge_indels(&p, h.nc);
    assert_eq!(nc.len(), 1);
    assert_eq!(nc[0].range, (20, 28));
    assert_eq!(nc[0].kind, InDelKind::Deletion);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_indel_no_interaction_hoisted() -> Result<(), Report> {
    // A parent indel disjoint from the child's indels is hoisted cleanly to u->N, leaving
    // N->v free of indels and N->c carrying only the child's own indel.
    let mut graph: GraphAncestral = nwk_read_str(NWK)?;
    let uv = find_edge_key(&graph, "U", "V").unwrap();
    let va = find_edge_key(&graph, "V", "A").unwrap();

    let partition = make_partition(
      &graph,
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
      let mut p = partition.write_arc();
      p.edges.get_mut(&uv).unwrap().indels = vec![parent_del.clone()];
      p.edges.get_mut(&va).unwrap().indels = vec![child_del.clone()];
    }
    let sparse = vec![partition];

    hoist_reverting_child(&mut graph, &sparse, &no_dense(), uv, va)?;

    let h = Hoisted::locate(&graph, "V", "A");
    let p = sparse[0].read_arc();
    assert_eq!(edge_indels(&p, h.un), vec![parent_del]);
    assert_eq!(edge_indels(&p, h.nv), Vec::<InDel>::new());
    assert_eq!(edge_indels(&p, h.nc), vec![child_del]);
    Ok(())
  }

  // Bifurcating root: root -> {V, S}. V groups two children; S is the sibling leaf. A site
  // where V equals the root but S differs carries its substitution on the sibling edge root->S,
  // and a child of V that makes the same change reverts it across the root.
  const NWK_BIFURCATING: &str = "((C1:0.1,C2:0.1)V:0.1,S:0.1)root:0.0;";

  #[test]
  fn test_hoist_reversions_slide_moves_sibling_sub_to_parent() -> Result<(), Report> {
    // Sibling edge root->S carries A3G (root=A, S=G); V's parent edge is empty at 3 (V=A);
    // child V->C1 carries A3G (the same change). The slide re-roots the site onto the
    // sibling's state: root becomes G, the sibling edge empties, and the parent edge gains the
    // inverse G3A. The root's clamped MAP sequence moves with root_sequence.
    let graph: GraphAncestral = nwk_read_str(NWK_BIFURCATING)?;
    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let root_v = find_edge_key(&graph, "root", "V").unwrap();
    let root_s = find_edge_key(&graph, "root", "S").unwrap();
    let v_c1 = find_edge_key(&graph, "V", "C1").unwrap();

    let partition = make_partition(
      &graph,
      0,
      100,
      &[
        ("root", "S", vec![sub(b'A', 3, b'G')]),
        ("V", "C1", vec![sub(b'A', 3, b'G')]),
      ],
    );
    let sparse = vec![partition];

    slide_bifurcating_root_for_child(&sparse, root_key, root_v, root_s, v_c1)?;

    let p = sparse[0].read_arc();
    assert_eq!(p.root_sequence[3], c(b'G'));
    assert_eq!(p.nodes[&root_key].seq.sequence[3], c(b'G'));
    assert_eq!(edge_subs(&p, root_s), Vec::<Sub>::new());
    assert_eq!(edge_subs(&p, root_v), vec![sub(b'G', 3, b'A')]);
    Ok(())
  }

  #[test]
  fn test_hoist_reversions_slide_then_hoist_removes_cross_root_reversion() -> Result<(), Report> {
    // The slide exposes the reversion on V's parent edge; the standard hoist then removes it.
    // Two mutations before (sibling A3G plus the child's A3G), one after (delta = -1): a single
    // G3A on the branch to the A-state V, exactly the reroot-invariant parsimony cost of the
    // site. The slide alone is count-neutral, so the reduction comes from the hoist it enables.
    let mut graph: GraphAncestral = nwk_read_str(NWK_BIFURCATING)?;
    let root_key = find_node_key_by_name(&graph, "root").unwrap();
    let root_v = find_edge_key(&graph, "root", "V").unwrap();
    let root_s = find_edge_key(&graph, "root", "S").unwrap();
    let v_c1 = find_edge_key(&graph, "V", "C1").unwrap();

    let partition = make_partition(
      &graph,
      0,
      100,
      &[
        ("root", "S", vec![sub(b'A', 3, b'G')]),
        ("V", "C1", vec![sub(b'A', 3, b'G')]),
      ],
    );
    let sparse = vec![partition];

    let before = helpers::total_subs(&graph, &sparse[0].read_arc());
    slide_bifurcating_root_for_child(&sparse, root_key, root_v, root_s, v_c1)?;
    let after_slide = helpers::total_subs(&graph, &sparse[0].read_arc());
    hoist_reverting_child(&mut graph, &sparse, &no_dense(), root_v, v_c1)?;
    let after_hoist = helpers::total_subs(&graph, &sparse[0].read_arc());

    assert_eq!(before, 2);
    assert_eq!(after_slide, 2); // the slide moves a substitution, it does not remove one
    assert_eq!(after_hoist, 1);
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition};
    use maplit::btreemap;
    use treetime_primitives::{AsciiChar, Seq, seq};

    pub fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(c(reff), pos, c(qry)).unwrap()
    }

    pub fn no_dense() -> Vec<Arc<RwLock<PartitionMarginalDense>>> {
      vec![]
    }

    /// Edge keys of the three edges the hoist produces, located by the child/node names.
    pub struct Hoisted {
      pub v: GraphNodeKey,
      pub un: GraphEdgeKey,
      pub nv: GraphEdgeKey,
      pub nc: GraphEdgeKey,
    }

    impl Hoisted {
      pub fn locate(graph: &GraphAncestral, v_name: &str, c_name: &str) -> Self {
        let v = find_node_key_by_name(graph, v_name).unwrap();
        let c = find_node_key_by_name(graph, c_name).unwrap();
        let nv = single_inbound(graph, v);
        let n = graph.get_source_node_key(nv).unwrap();
        let un = single_inbound(graph, n);
        let nc = single_inbound(graph, c);
        Self { v, un, nv, nc }
      }
    }

    fn single_inbound(graph: &GraphAncestral, node_key: GraphNodeKey) -> GraphEdgeKey {
      let node = graph.get_node(node_key).unwrap();
      let node = node.read_arc();
      match node.inbound() {
        [edge_key] => *edge_key,
        other => panic!("expected exactly one inbound edge, found {}", other.len()),
      }
    }

    pub fn edge_subs(partition: &PartitionMarginalSparse, edge_key: GraphEdgeKey) -> Vec<Sub> {
      partition.edges[&edge_key].fitch_subs().to_vec()
    }

    pub fn edge_indels(partition: &PartitionMarginalSparse, edge_key: GraphEdgeKey) -> Vec<InDel> {
      partition.edges[&edge_key].indels.clone()
    }

    pub fn total_subs(graph: &GraphAncestral, partition: &PartitionMarginalSparse) -> usize {
      graph
        .get_edges()
        .iter()
        .filter_map(|e| partition.edges.get(&e.read_arc().key()))
        .map(|e| e.fitch_subs().len())
        .sum()
    }

    pub fn branch_length(graph: &GraphAncestral, edge_key: GraphEdgeKey) -> f64 {
      graph
        .get_edge(edge_key)
        .and_then(|e| e.read_arc().payload().read_arc().branch_length())
        .unwrap()
    }

    pub fn make_partition(
      graph: &GraphAncestral,
      index: usize,
      length: usize,
      edge_mutations: &[(&str, &str, Vec<Sub>)],
    ) -> Arc<RwLock<PartitionMarginalSparse>> {
      let mut partition = PartitionMarginalSparse {
        index,
        gtr: jc69(JC69Params::default()).unwrap(),
        alphabet: Alphabet::new(AlphabetName::Nuc).unwrap(),
        length,
        nodes: btreemap! {},
        edges: btreemap! {},
        root_sequence: seq![],
      };

      let mut ref_seq: Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
      for (_, _, subs) in edge_mutations {
        for s in subs {
          if s.pos() < length {
            ref_seq[s.pos()] = s.reff();
          }
        }
      }
      partition.root_sequence = ref_seq.clone();

      for node in graph.get_nodes() {
        let key = node.read_arc().key();
        let mut node_part = SparseNodePartition::empty(&partition.alphabet);
        node_part.seq.sequence = ref_seq.clone();
        partition.nodes.insert(key, node_part);
      }

      for (source, target, subs) in edge_mutations {
        let edge_key =
          find_edge_key(graph, source, target).unwrap_or_else(|| panic!("edge {source}->{target} missing"));
        partition
          .edges
          .insert(edge_key, SparseEdgePartition::with_fitch_subs(subs.clone()));
      }

      Arc::new(RwLock::new(partition))
    }
  }
}
