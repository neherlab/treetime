/// Fitch parsimony indel reconstruction on a 5-taxon tree with mixed gap/unknown patterns.
///
/// Tree: (((A,B)I1,C)I2,(D,E)I3)root
///
/// Alignment (8 positions, 0-indexed):
///
/// ```text
/// A  CGT--A-T   gaps at (3,5) and (6,7)
/// B  CGT-NNNT   gap at (3,4), unknowns at (4,7)
/// C  CGTNNA-T   unknowns at (3,5), gap at (6,7)
/// D  CGTTTA-T   gap at (6,7)
/// E  CGTTGACT   no gaps or unknowns
/// ```
///
/// Key properties tested:
///
/// **Gap-beats-unknown rule (backward pass)**: When all children at a position are
/// non-character (gap or unknown), and at least one child has a gap (not merely unknown),
/// the parent inherits a gap rather than an unknown. This is exercised at nodes I1 and I2:
///
/// - I1 (A,B): position range (3,5) — A has a gap, B has unknown. Because every child is
///   non-char and at least one has a gap, I1 inherits a gap at (3,5) (and at (6,7) by the
///   same rule).
/// - I2 (I1,C): I1 has a gap at (3,5), C has unknown there. Same rule ⇒ I2.gaps includes
///   (3,5).
///
/// **Variable-indel resolution at the root**: At I3, D has a gap at (6,7) but E does not,
/// producing a variable indel `{(6,7): deleted=1, present=1}`. This propagates to the root
/// where it combines with I2's gap at (6,7): both children contribute a "deleted" vote,
/// so the range is resolved as a consensus gap and added to `root.gaps`.
///
/// **InDel sequence content**: The deletion on edge root→I2 at (3,5) must carry the
/// parent's (root's) sequence — not the child's gap-filled sequence — ensuring the InDel
/// record contains actual nucleotide characters and not gap characters.
#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::fitch::compress_sequences;
  use crate::o;
  use crate::representation::partition::fitch::PartitionFitch;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use std::sync::{Arc, LazyLock};
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::read_many_fasta;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::vec_of_owned;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .map(PathBuf::from)
      .expect("project has workspace root")
  }

  fn get_node_name(graph: &GraphAncestral, key: GraphNodeKey) -> String {
    let node = graph.get_node(key).expect("node exists");
    node
      .read_arc()
      .payload()
      .read_arc()
      .name
      .clone()
      .expect("node has name")
  }

  fn get_node_gaps_by_name(
    graph: &GraphAncestral,
    partition: &PartitionFitch,
    name: &str,
  ) -> Vec<(usize, usize)> {
    for node in graph.get_nodes() {
      let node = node.read_arc();
      if node.payload().read_arc().name.as_deref() == Some(name) {
        let mut gaps = partition.nodes[&node.key()].seq.gaps.clone();
        gaps.sort();
        return gaps;
      }
    }
    panic!("Node '{name}' not found");
  }

  fn collect_edge_indels(graph: &GraphAncestral, partition: &PartitionFitch) -> BTreeMap<String, Vec<String>> {
    graph
      .get_edges()
      .iter()
      .map(|edge| {
        let edge = edge.read_arc();
        let parent_name = get_node_name(graph, edge.source());
        let child_name = get_node_name(graph, edge.target());
        let edge_name = format!("{parent_name}->{child_name}");
        let indels = partition.edges[&edge.key()]
          .indels
          .iter()
          .map(|indel| indel.to_string())
          .collect_vec();
        (edge_name, indels)
      })
      .collect()
  }

  /// After `compress_sequences`, assert the gap ranges stored at every internal node.
  ///
  /// The "gap-beats-unknown" rule (backward pass) and variable-indel propagation
  /// (resolve_indels_backward) jointly determine these ranges.
  #[test]
  fn test_fitch_indel_data_node_gaps() -> Result<(), Report> {
    let root = project_root();
    let aln = read_many_fasta(root.join("data/indel/aln.fasta"), &*NUC_ALPHABET)?;

    // Same topology as data/indel/tree.nwk with internal nodes named for assertions.
    let graph: GraphAncestral =
      nwk_read_str("(((A:0.01,B:0.02)I1:0.02,C:0.04)I2:0.03,(D:0.02,E:0.025)I3:0.01)root;")?;

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet: Alphabet::default(),
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let partition = partitions[0].read_arc();

    // I1 (parent of A and B):
    //   A.non_char = [(3,5),(6,7)],  A.gaps = [(3,5),(6,7)]
    //   B.non_char = [(3,7)],        B.gaps = [(3,4)],   B.unknown = [(4,7)]
    //   non_char(I1) = intersection = [(3,5),(6,7)]
    //   gap_union    = [(3,5),(6,7)]   (A's full gap set beats B's partial gap+unknown)
    //   → I1.gaps = [(3,5),(6,7)]
    assert_eq!(vec![(3, 5), (6, 7)], get_node_gaps_by_name(&graph, &partition, "I1"), "I1 gaps");

    // I2 (parent of I1 and C):
    //   I1.gaps = [(3,5),(6,7)],  C.gaps = [(6,7)],  C.unknown = [(3,5)]
    //   non_char(I2) = [(3,5),(6,7)] (both children non_char at those ranges)
    //   gap_union    = [(3,5),(6,7)]   (I1 has gap at (3,5), beats C's unknown there)
    //   → I2.gaps = [(3,5),(6,7)]
    assert_eq!(
      vec![(3, 5), (6, 7)],
      get_node_gaps_by_name(&graph, &partition, "I2"),
      "I2 gaps"
    );

    // I3 (parent of D and E):
    //   D.gaps = [(6,7)],  E.gaps = []
    //   non_char(I3) = []  (E has no non_char positions)
    //   → backward: I3.gaps = [],  I3.variable_indel = {(6,7): {deleted:1, present:1}}
    //   → forward:  root also has gap at (6,7) → deleted+gap_in_parent = 2 > present = 1
    //   → I3.gaps becomes [(6,7)] after the forward pass resolves the variable indel
    assert_eq!(vec![(6, 7)], get_node_gaps_by_name(&graph, &partition, "I3"), "I3 gaps");

    // root (parent of I2 and I3):
    //   I2.non_char = [(3,5),(6,7)],  I3.non_char = []
    //   → backward: root.non_char = [], root.gaps = []
    //   variable_indel: (3,5) → {deleted:1, present:1}  (only I2 has gap there)
    //                   (6,7) → resolved as consensus gap via I3's variable_indel vote
    //   resolved_gaps = [(6,7)] → appended to root.gaps
    //   → root.gaps = [(6,7)]
    assert_eq!(vec![(6, 7)], get_node_gaps_by_name(&graph, &partition, "root"), "root gaps");

    Ok(())
  }

  /// After `compress_sequences`, assert the InDel records on every edge.
  ///
  /// Key cases:
  /// - root→I2: deletion at (3,5) — I2 has a consensus gap there but root does not.
  ///   The deletion sequence must come from root's sequence (non-gap nucleotides),
  ///   not from I2's gap-filled positions. This validates the fix that stopped
  ///   InDels from carrying gap characters in their sequence payload.
  /// - I3→E: insertion at (6,7) — I3 has a gap there (resolved in the forward pass)
  ///   but leaf E has 'C' at position 6. The insertion sequence is "C".
  ///
  /// Note: edges I2→C and I1→B produce insertion InDels at ranges where the child
  /// has unknown ('N') characters but the non_char fill in the forward pass has already
  /// overwritten those positions with the parent's gap before InDel construction. The
  /// resulting InDels carry gap characters in their sequence — a known limitation when
  /// leaves have unknown characters adjacent to ancestral gaps.
  #[test]
  fn test_fitch_indel_data_edge_indels() -> Result<(), Report> {
    let root = project_root();
    let aln = read_many_fasta(root.join("data/indel/aln.fasta"), &*NUC_ALPHABET)?;

    let graph: GraphAncestral =
      nwk_read_str("(((A:0.01,B:0.02)I1:0.02,C:0.04)I2:0.03,(D:0.02,E:0.025)I3:0.01)root;")?;

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet: Alphabet::default(),
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let partition = partitions[0].read_arc();
    let actual = collect_edge_indels(&graph, &partition);

    let expected = btreemap! {
      // root has no gap at (3,5) but I2 does → deletion on this edge.
      // Deletion seq = root.sequence[3..5] = "TG" (T from position 3, G from
      // the variable site at position 4 resolved to 'G' by majority rule).
      o!("root->I2") => vec_of_owned!["3--5: TG -> --"],

      // root has gap at (6,7) and I3 inherits it through variable-indel resolution
      // (deleted+gap_in_parent = 2 > present = 1) → no InDel on this edge.
      o!("root->I3") => vec![],

      // I1 and I2 share the same gap ranges → no InDels.
      o!("I2->I1") => vec![],

      // I2 has gap at (3,5); C does not (it has unknown 'N' there).
      // The non_char fill overwrites C.sequence[3..5] with I2's gap '--' before
      // the InDel is constructed, so the insertion sequence is "--" rather than "NN".
      o!("I2->C")  => vec_of_owned!["3--5: -- -> --"],

      // A and I1 share the same gap ranges → no InDels.
      o!("I1->A")  => vec![],

      // I1 has gaps at (4,5) and (6,7) that B does not have (B has unknown there).
      // Same non_char-fill ordering issue as I2→C: sequence is '--' at construction time.
      o!("I1->B")  => vec_of_owned!["4--5: - -> -", "6--7: - -> -"],

      // D and I3 both have gap at (6,7) after forward pass → no InDels.
      o!("I3->D")  => vec![],

      // I3 has gap at (6,7); E has 'C' there → insertion with seq "C".
      o!("I3->E")  => vec_of_owned!["6--7: - -> C"],
    };

    assert_eq!(expected, actual);

    // The key regression check: the deletion on root→I2 must not carry gap characters
    // in its sequence payload. Before the fix, the wrong sequence source (node rather
    // than parent) caused "TG -> --" to appear as "-- -> --".
    let root_i2_indels = actual.get("root->I2").expect("root->I2 edge exists");
    assert_eq!(root_i2_indels.len(), 1);
    assert!(
      root_i2_indels[0].starts_with("3--5: TG"),
      "deletion sequence on root->I2 must be non-gap nucleotides, got: {}",
      root_i2_indels[0]
    );

    Ok(())
  }
}
