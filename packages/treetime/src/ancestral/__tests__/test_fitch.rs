#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::{
    ancestral_reconstruction_fitch, attach_seqs_to_graph, compress_sequences, fitch_backward, fitch_forward,
  };
  use crate::ancestral::marginal::update_marginal;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::o;
  use crate::partition::fitch::PartitionFitch;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::traits::{PartitionBranchOps, PartitionRerootOps};
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use crate::seq::composition::Composition;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::{Arc, LazyLock};
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::reroot::{RerootChanges, apply_reroot_topology, remove_node_if_trivial, split_edge};
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AsciiChar;
  use treetime_utils::io::json::{JsonPretty, json_write_str};
  use treetime_utils::sync::mutex::unwrap_arc_rwlock;
  use treetime_utils::vec_of_owned;

  /// Retrieve the name of a graph node by its key. Panics if the node is missing or unnamed.
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

  /// Collect substitution mutations on every edge, keyed by "parent->child" label.
  ///
  /// Each edge's substitutions are formatted as strings (e.g. "A1G" for A->G at position 1).
  /// Used to verify that compression (backward + forward + cleanup) correctly placed
  /// mutations on branches.
  fn collect_edge_subs(graph: &GraphAncestral, partition: &PartitionFitch) -> BTreeMap<String, Vec<String>> {
    graph
      .get_edges()
      .iter()
      .map(|edge| {
        let edge = edge.read_arc();
        let parent_name = get_node_name(graph, edge.source());
        let child_name = get_node_name(graph, edge.target());
        let edge_name = format!("{parent_name}->{child_name}");
        let subs = partition.edges[&edge.key()]
          .fitch_subs()
          .iter()
          .map(|sub| sub.to_string())
          .collect_vec();
        (edge_name, subs)
      })
      .collect()
  }

  /// Return the sorted list of alignment positions that are variable at the root node.
  ///
  /// After Fitch's backward pass, variable positions are those where the intersection of
  /// child state sets was either empty (union taken) or contained more than one state
  /// (ambiguous intersection). Positions where the intersection was a singleton are resolved
  /// immediately and do not appear as variable.
  fn get_root_variable_positions(graph: &GraphAncestral, partition: &PartitionFitch) -> Vec<usize> {
    let root = graph.get_exactly_one_root().expect("graph has exactly one root");
    let root_key = root.read_arc().key();
    partition.nodes[&root_key]
      .seq
      .fitch
      .variable
      .keys()
      .copied()
      .collect_vec()
  }

  /// Return the Fitch state sets at each variable position of the root node.
  ///
  /// After the backward pass, each variable position holds the set of candidate nucleotide
  /// states. The set is either a multi-element intersection (children partially agree) or the
  /// union of child state sets (children share no common state, implying a state change).
  /// Singleton intersections (all children agree on one state) are resolved immediately and
  /// do not appear in this map.
  fn get_root_state_sets(graph: &GraphAncestral, partition: &PartitionFitch) -> BTreeMap<usize, String> {
    let root = graph.get_exactly_one_root().expect("graph has exactly one root");
    let root_key = root.read_arc().key();
    partition.nodes[&root_key]
      .seq
      .fitch
      .variable
      .iter()
      .map(|(&pos, states)| (pos, states.to_string()))
      .collect()
  }

  /// Look up a node by name and return its variable positions after the Fitch backward pass.
  /// Panics if no node with the given name exists.
  fn get_node_variable_positions_by_name(graph: &GraphAncestral, partition: &PartitionFitch, name: &str) -> Vec<usize> {
    for node in graph.get_nodes() {
      let node = node.read_arc();
      let node_name = node.payload().read_arc().name.clone();
      if node_name.as_deref() == Some(name) {
        return partition.nodes[&node.key()]
          .seq
          .fitch
          .variable
          .keys()
          .copied()
          .collect_vec();
      }
    }
    panic!("Node {name} not found");
  }

  /// Return the reconstructed sequence at the root node as a string.
  ///
  /// After the backward pass, unresolved variable positions are marked with '~'.
  /// After the forward pass, all positions are resolved to concrete nucleotides.
  fn get_root_seq(graph: &GraphAncestral, partition: &PartitionFitch) -> String {
    let root = graph.get_exactly_one_root().expect("graph has exactly one root");
    let root_key = root.read_arc().key();
    partition.nodes[&root_key].seq.sequence.as_str().to_owned()
  }

  fn get_internal_sequences(graph: &GraphAncestral, partition: &PartitionFitch) -> BTreeMap<String, String> {
    graph
      .get_internal_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let node_name = node.payload().read_arc().name.clone().unwrap();
        let sequence = partition.nodes[&node.key()].seq.sequence.as_str().to_owned();
        (node_name, sequence)
      })
      .collect()
  }

  /// Collect insertion/deletion (indel) mutations on every edge, keyed by "parent->child" label.
  ///
  /// Each indel is formatted as a range with parent and child states (e.g. "12--13: T -> -"
  /// for a deletion of T at positions 12-13). Used to verify that the Fitch algorithm correctly
  /// identifies gap openings, extensions, and insertions relative to the ancestral sequence.
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

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  /// Fitch maximum parsimony reconstruction on a balanced 4-taxon binary tree (Fitch 1971).
  ///
  /// Tree topology: ((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01
  ///
  /// Verifies that the full Fitch pipeline (backward + forward + cleanup, then reconstruction
  /// traversal) produces the correct ancestral sequences at internal nodes (root, AB, CD).
  /// The alignment contains substitutions, ambiguous bases (N, R), and gap patterns to
  /// exercise all code paths.
  ///
  /// Only internal node sequences are emitted (leaves excluded via `include_leaves=false`).
  /// Expected sequences are hand-derived from the Fitch intersection/union rule applied
  /// bottom-up, then resolved top-down: at the root, the alphabetically first state from
  /// the candidate set is chosen (no parent to prefer); at non-root nodes, the parent
  /// state is preferred when present in the child's state set.
  #[test]
  fn test_ancestral_reconstruction_fitch() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
        >root
        ACAGCCATGTATTG--
        >AB
        ACATCCCTGTA--G--
        >CD
        CCGGCCATGTATTG--
      "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions_parsimony = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions_parsimony, &aln)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_fitch(&graph, false, &partitions_parsimony, |node, seq| {
      actual.insert(node.payload.name.clone(), seq.to_string());
      Ok(())
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  /// Same Fitch reconstruction as `test_ancestral_reconstruction_fitch`, but with
  /// `include_leaves=true` so that both internal and leaf node sequences are emitted.
  ///
  /// Verifies that leaf sequences pass through unchanged (the algorithm must not alter
  /// observed data) and that internal node sequences match the expected reconstruction.
  #[test]
  fn test_ancestral_reconstruction_fitch_with_leaves() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
        >root
        ACAGCCATGTATTG--
        >AB
        ACATCCCTGTA--G--
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >CD
        CCGGCCATGTATTG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions_parsimony = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions_parsimony, &aln)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_fitch(&graph, true, &partitions_parsimony, |node, seq| {
      actual.insert(node.payload.name.clone(), seq.to_string());
      Ok(())
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_compress_sequences_retains_internal_exact_sequences() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        RCGTACGT
        >B
        GCGTACGT
        >C
        GCGTACGT
        >D
        GCGTACGT
      "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let alphabet = Alphabet::default();

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let partition = partitions[0].read_arc();
    let actual = get_internal_sequences(&graph, &partition);
    let expected = btreemap! {
      o!("AB") => o!("GCGTACGT"),
      o!("CD") => o!("GCGTACGT"),
      o!("root") => o!("GCGTACGT"),
    };

    assert_eq!(expected, actual);
    Ok(())
  }

  /// Verify substitution and indel mutations placed on individual edges by Fitch's algorithm.
  ///
  /// Uses the same 4-taxon tree and alignment as `test_ancestral_reconstruction_fitch`.
  /// After running `compress_sequences()` (backward + forward + cleanup), inspects each
  /// edge for:
  /// - Substitutions: point mutations between parent and child states (e.g. "C6G")
  /// - Indels: gap openings/closings with position ranges (e.g. "12--13: T -> -")
  ///
  /// This tests the mutation-placement logic in the forward pass, verifying that the diff
  /// between parent and child sequences is correctly decomposed into substitution and
  /// indel events.
  #[test]
  fn test_fitch_internals() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let partition = partitions[0].read_arc();

    // Verify substitutions on edges
    let actual_subs = collect_edge_subs(&graph, &partition);
    let expected_subs = btreemap! {
      o!("AB->A")    => vec_of_owned!["C6G", "T8C"],
      o!("AB->B")    => vec_of_owned!["A1G"],
      o!("root->AB") => vec_of_owned!["G4T", "A7C"],
      o!("CD->C")    => vec_of_owned!["C6G"],
      o!("CD->D")    => vec_of_owned!["C1T", "A7G"],
      o!("root->CD") => vec_of_owned!["A1C", "A3G"],
    };
    assert_eq!(expected_subs, actual_subs);

    // Verify indels on edges
    let actual_indels = collect_edge_indels(&graph, &partition);
    let expected_indels = btreemap! {
      o!("AB->A")     => vec_of_owned!["14--16: -- -> AC"],
      o!("AB->B")     => vec![],
      o!("root->AB")  => vec_of_owned!["11--13: TT -> --"],
      o!("CD->C")     => vec![],
      o!("CD->D")     => vec![],
      o!("root->CD")  => vec![],
    };
    assert_eq!(expected_indels, actual_indels);

    Ok(())
  }

  /// Fitch algorithm with complex gap (indel) patterns on a 5-taxon asymmetric tree.
  ///
  /// Tree: ((A,B)AB,(C,(D,E)DE)CDE)root
  ///
  /// Standard Fitch parsimony (Fitch 1971) operates on single characters. This
  /// implementation extends it with range-based gap tracking for multi-position insertions
  /// and deletions. Tests three challenging indel scenarios via `compress_sequences()`
  /// (backward + forward + cleanup):
  /// - Overlapping deletions: leaves A and B have deletions at different but overlapping
  ///   positions, requiring the algorithm to reconstruct which gaps are ancestral vs derived.
  /// - Root-level deletion: the root itself carries a gap that propagates to subtrees.
  /// - Variable inserted sequence: positions 2-3 in subtree DE are inserted (absent in
  ///   the rest of the tree), and the inserted characters differ between D and E,
  ///   requiring both indel and substitution tracking within the insertion.
  #[test]
  fn test_fitch_complex_gaps() -> Result<(), Report> {
    // Test cases: a) deletions overlap, b) root has a deletion, c) inserted sequence is variable
    // In DE, position 3 is inserted, but it varies in D and E
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        NC--G
        >B
        T--AG
        >C
        TR-TG
        >D
        TGTTG
        >E
        TGCCG
      "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,(D:0.05,E:0.03)DE:0.01)CDE:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let partition = partitions[0].read_arc();

    // Verify substitutions on edges
    let actual_subs = collect_edge_subs(&graph, &partition);
    let expected_subs = btreemap! {
      o!("AB->A")     => vec![],
      o!("AB->B")     => vec![],
      o!("root->AB")  => vec![],
      o!("CDE->C")    => vec![],
      o!("CDE->DE")   => vec![],
      o!("DE->D")     => vec![o!("C3T")],
      o!("DE->E")     => vec![o!("T4C")],
      o!("root->CDE") => vec![o!("C2G"), o!("A4T")],
    };
    assert_eq!(expected_subs, actual_subs);

    // Verify indels on edges
    let actual_indels = collect_edge_indels(&graph, &partition);
    let expected_indels = btreemap! {
      o!("AB->A")     => vec![o!("3--4: A -> -")],
      o!("AB->B")     => vec![o!("1--2: C -> -")],
      o!("root->AB")  => vec![],
      o!("CDE->C")    => vec![],
      o!("CDE->DE")   => vec![o!("2--3: - -> C")],
      o!("DE->D")     => vec![],
      o!("DE->E")     => vec![],
      o!("root->CDE") => vec![],
    };
    assert_eq!(expected_indels, actual_indels);

    Ok(())
  }

  /// Fitch algorithm on a tree with a polytomy (multifurcation).
  ///
  /// Tree: ((A,B)AB,(C,D,E)CDE)root - node CDE has 3 children instead of 2.
  ///
  /// Fitch's algorithm generalizes from binary to n-ary nodes: the backward pass computes
  /// the intersection of all n child state sets. If the intersection is empty, it takes
  /// the union (implying at least one state change). The forward pass resolves ambiguities
  /// using the parent state when it is present in the child's state set.
  ///
  /// Uses the same alignment as `test_fitch_complex_gaps` but with a different tree topology
  /// (D and E are direct children of CDE rather than grouped under a DE intermediate node).
  /// This changes which substitutions and indels appear on which edges compared to the
  /// binary tree version.
  #[test]
  fn test_fitch_polytomy() -> Result<(), Report> {
    // Test polytomy (node with more than 2 children)
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        NC--G
        >B
        T--AG
        >C
        TR-TG
        >D
        TGTTG
        >E
        TGCCG
      "#},
      &*NUC_ALPHABET,
    )?;

    // CDE is a polytomy with 3 children: C, D, E
    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.05,E:0.03)CDE:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let partition = partitions[0].read_arc();

    // Verify substitutions on edges
    let actual_subs = collect_edge_subs(&graph, &partition);
    let expected_subs = btreemap! {
      o!("AB->A")     => vec![],
      o!("AB->B")     => vec![],
      o!("root->AB")  => vec![],
      o!("CDE->C")    => vec![o!("C4T")],
      o!("CDE->D")    => vec![o!("C4T")],
      o!("CDE->E")    => vec![],
      o!("root->CDE") => vec![o!("C2G"), o!("A4C")],
    };
    assert_eq!(expected_subs, actual_subs);

    // Verify indels on edges
    let actual_indels = collect_edge_indels(&graph, &partition);
    let expected_indels = btreemap! {
      o!("AB->A")     => vec![o!("3--4: A -> -")],
      o!("AB->B")     => vec![o!("1--2: C -> -")],
      o!("root->AB")  => vec![],
      o!("CDE->C")    => vec![],
      o!("CDE->D")    => vec![o!("2--3: - -> T")],
      o!("CDE->E")    => vec![o!("2--3: - -> C")],
      o!("root->CDE") => vec![],
    };
    assert_eq!(expected_indels, actual_indels);

    Ok(())
  }

  /// Inspect intermediate state of Fitch's algorithm between the backward and forward passes.
  ///
  /// Fitch's algorithm (Fitch 1971) has two distinct phases:
  ///
  /// **Backward pass** (post-order, leaves to root): At each internal node, compute the
  /// intersection of child state sets. Three outcomes per position:
  /// - Singleton intersection (all children agree on one state): resolved immediately.
  /// - Multi-element intersection (children partially agree): stored as variable, marked '~'.
  /// - Empty intersection (children share no state, implying a state change): the union is
  ///   stored as variable, marked '~'.
  ///
  /// **Forward pass** (pre-order, root to leaves): Resolve variable state sets to concrete
  /// nucleotides. At the root (no parent), the alphabetically first state from the candidate
  /// set is chosen via `StateSet::get_one()`. At non-root nodes, the parent state is chosen
  /// when present in the child's state set. After this pass, all '~' markers are replaced.
  ///
  /// This test runs each pass separately and verifies:
  /// - After backward: correct variable positions identified, state sets contain the
  ///   expected candidate nucleotides, root sequence has '~' at unresolved positions.
  /// - After forward: same variable positions are tracked, but root sequence is fully resolved.
  #[test]
  fn test_fitch_backward_state() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    // Run backward pass only
    attach_seqs_to_graph(&graph, &partitions, &aln)?;
    fitch_backward(&graph, &partitions);

    {
      let partition = partitions[0].read_arc();

      // After backward: variable positions identified, sequence has '~' markers
      let variable_positions = get_root_variable_positions(&graph, &partition);
      assert_eq!(vec![0, 2, 3, 5, 6], variable_positions);

      let root_seq = get_root_seq(&graph, &partition);
      assert_eq!("~C~~C~~TGTATTGAC", root_seq);

      // Verify state sets at each variable position (which characters are possible)
      let state_sets = get_root_state_sets(&graph, &partition);
      assert_eq!(
        btreemap! {
          0 => o!("{A, C, G, T}"),  // union: AB={A,G}, CD={C,T}, intersection empty
          2 => o!("{A, G}"),        // union: AB=A (resolved), CD=G (resolved), intersection empty
          3 => o!("{G, T}"),        // union: AB=T (resolved), CD=G (resolved), intersection empty
          5 => o!("{C, G}"),        // intersection: AB={C,G}, CD={C,G}, intersection={C,G} (ambiguous)
          6 => o!("{A, C, G}"),     // union: AB=C (resolved), CD={A,G}, intersection empty
        },
        state_sets
      );

      // Verify internal node AB has variable positions (specific values depend on leaf data)
      let ab_vars = get_node_variable_positions_by_name(&graph, &partition, "AB");
      assert!(!ab_vars.is_empty(), "AB should have variable positions");

      // Verify internal node CD has variable positions
      let cd_vars = get_node_variable_positions_by_name(&graph, &partition, "CD");
      assert!(!cd_vars.is_empty(), "CD should have variable positions");
    }

    // Run forward pass
    fitch_forward(&graph, &partitions)?;

    {
      let partition = partitions[0].read_arc();

      // After forward: variable positions still tracked, sequence resolved
      let variable_positions = get_root_variable_positions(&graph, &partition);
      assert_eq!(vec![0, 2, 3, 5, 6], variable_positions);

      let root_seq = get_root_seq(&graph, &partition);
      assert_eq!("ACAGCCATGTATTG--", root_seq);
    }

    Ok(())
  }

  /// Verify that the sparse partition remains consistent after rerooting on branch AB->A.
  ///
  /// Uses the same 4-taxon tree and alignment as `test_fitch_internals`. After Fitch
  /// reconstruction, converts to a sparse partition, then places the new root at the
  /// midpoint of the AB->A edge using `split_edge`. The old root is kept as a trivial node
  /// (no merge). Verifies:
  ///
  /// - The root_sequence after reroot equals the AB node's sequence (new root sits on
  ///   the AB->A branch, so it inherits AB's ancestral state with no further mutations).
  /// - The child-side edge (new_root->A) retains the original AB->A substitutions and
  ///   indels unchanged (mutations still describe AB-state to A-state).
  /// - The parent-side edge (new_root->AB) is empty (no mutations between the split point
  ///   and AB, since the split was placed at 0.0 relative to AB).
  /// - The inverted root->AB edge (now AB->root) has its substitutions and indels inverted:
  ///   reff and qry swapped, deletion flag toggled.
  #[test]
  fn test_fitch_reroot_sparse_on_branch_ab_to_a() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let alphabet = Alphabet::default();

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let fitch = unwrap_arc_rwlock(partitions.into_iter().next().unwrap())?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let mut sparse = fitch.into_marginal_sparse(gtr, &graph)?;

    // Find relevant node keys and the AB->A edge key
    let old_root_key = graph.get_exactly_one_root()?.read_arc().key();
    let ab_key = find_node_key_by_name(&graph, "AB").expect("AB node not found");
    let a_key = find_node_key_by_name(&graph, "A").expect("A node not found");

    let edge_ab_a_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == ab_key && e.target() == a_key
      })
      .map(|e| e.read_arc().key())
      .expect("AB->A edge not found");

    // Record original AB->A subs and indels before the split
    let orig_subs: Vec<String> = sparse.edges[&edge_ab_a_key]
      .fitch_subs()
      .iter()
      .map(|s| s.to_string())
      .collect();
    let orig_indels: Vec<String> = sparse.edges[&edge_ab_a_key]
      .indels
      .iter()
      .map(|i| i.to_string())
      .collect();

    // Find the root->AB edge key (for verifying inversion)
    let edge_root_ab_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == old_root_key && e.target() == ab_key
      })
      .map(|e| e.read_arc().key())
      .expect("root->AB edge not found");

    let orig_root_ab_subs: Vec<String> = sparse.edges[&edge_root_ab_key]
      .fitch_subs()
      .iter()
      .map(|s| s.to_string())
      .collect();
    let orig_root_ab_indels: Vec<String> = sparse.edges[&edge_root_ab_key]
      .indels
      .iter()
      .map(|i| i.to_string())
      .collect();

    // Split AB->A at 0.5 to insert new root node
    let split_info = split_edge(&mut graph, edge_ab_a_key, 0.5)?;
    let new_root_key = split_info.new_node_key;
    let parent_side_key = split_info.parent_side_edge_key; // AB -> new_root
    let child_side_key = split_info.child_side_edge_key; // new_root -> A

    // Invert graph topology: edges on path new_root -> old_root are inverted
    let inverted_edge_keys = apply_reroot_topology(&mut graph, old_root_key, new_root_key)?;

    // Keep old root as trivial node (no merge)
    let changes = RerootChanges {
      edge_split: Some(split_info),
      edge_merge: None,
      inverted_edge_keys,
    };

    sparse.apply_reroot(&changes)?;

    // --- Verify root_sequence ---
    // New root sits on AB->A (closer to AB side, split at 0.5 with empty parent-side).
    // The only mutations on the path from old_root to new_root are via root->AB (now AB->root).
    // After inversion, root_seq should reflect AB's state: ACATCCCTGTA--G--
    let expected_root_seq = "ACATCCCTGTA--G--";
    let actual_root_seq = sparse.root_sequence.as_str();
    assert_eq!(
      expected_root_seq, actual_root_seq,
      "root_sequence after reroot should equal AB's sequence"
    );

    // --- Verify child-side edge (new_root->A) ---
    // Should carry the original AB->A subs and indels unchanged
    let child_edge = &sparse.edges[&child_side_key];
    let child_subs: Vec<String> = child_edge.fitch_subs().iter().map(|s| s.to_string()).collect();
    let child_indels: Vec<String> = child_edge.indels.iter().map(|i| i.to_string()).collect();
    assert_eq!(
      orig_subs, child_subs,
      "child-side edge subs should equal original AB->A subs"
    );
    assert_eq!(
      orig_indels, child_indels,
      "child-side edge indels should equal original AB->A indels"
    );

    // --- Verify parent-side edge (new_root->AB, was AB->new_root before inversion) ---
    // Should be empty (no mutations between split point and AB)
    let parent_edge = &sparse.edges[&parent_side_key];
    assert!(
      parent_edge.fitch_subs().is_empty(),
      "parent-side edge should have no subs"
    );
    assert!(parent_edge.indels.is_empty(), "parent-side edge should have no indels");

    // --- Verify inverted root->AB edge (now AB->root) ---
    // Original root->AB: subs=[G4T, A7C], indels=[11--13: TT->--]
    // After inversion: subs=[T4G, C7A], indels=[11--13: --->TT]
    let inv_edge = &sparse.edges[&edge_root_ab_key];
    let inv_subs: Vec<String> = inv_edge.fitch_subs().iter().map(|s| s.to_string()).collect();
    let inv_indels: Vec<String> = inv_edge.indels.iter().map(|i| i.to_string()).collect();
    assert_eq!(
      vec_of_owned!["T4G", "C7A"],
      inv_subs,
      "inverted AB->root subs should have reff/qry swapped"
    );
    assert_eq!(
      vec_of_owned!["11--13: -- -> TT"],
      inv_indels,
      "inverted AB->root indel should be toggled to insertion"
    );

    // Sanity: original root->AB subs were [G4T, A7C]
    assert_eq!(vec_of_owned!["G4T", "A7C"], orig_root_ab_subs);
    assert_eq!(vec_of_owned!["11--13: TT -> --"], orig_root_ab_indels);

    // --- G1: new root node composition matches root_sequence ---
    // ACATCCCTGTA--G--: A=3 C=4 G=2 T=3 -=4
    let root_node = &sparse.nodes[&new_root_key];
    let c = AsciiChar::from_byte_unchecked;
    #[rustfmt::skip]
    let expected_comp = Composition::from_counts(
      btreemap! {
        c(b'-') => 4, c(b'A') => 3, c(b'B') => 0, c(b'C') => 4,
        c(b'D') => 0, c(b'G') => 2, c(b'H') => 0, c(b'K') => 0,
        c(b'M') => 0, c(b'N') => 0, c(b'R') => 0, c(b'S') => 0,
        c(b'T') => 3, c(b'V') => 0, c(b'W') => 0, c(b'Y') => 0,
      },
      c(b'-'),
    );
    assert_eq!(
      expected_comp, root_node.seq.composition,
      "new root node composition should match root_sequence character counts"
    );

    // --- G2: gap and non_char ranges match root_sequence ---
    // Root sequence ACATCCCTGTA--G--: gaps at positions 11-12 and 14-15
    assert_eq!(
      vec![(11, 13), (14, 16)],
      root_node.seq.gaps,
      "new root gaps should cover gap positions in root_sequence"
    );
    assert!(
      root_node.seq.unknown.is_empty(),
      "new root should have no unknown (N) positions"
    );
    assert_eq!(
      root_node.seq.gaps, root_node.seq.non_char,
      "non_char should equal gaps when there are no N positions"
    );

    // --- G3: edge_effective_length excludes non_char positions ---
    // Root non_char: positions 11-12, 14-15 (4 positions)
    // Child A non_char: gaps at 11-12 + N at 8-9 (from original sequence ACATCGCCNNA--GAC)
    // Union: positions 8-9, 11-12, 14-15 = 6 positions
    // Effective length = 16 - 6 = 10
    let effective = sparse.edge_effective_length(&graph, child_side_key)?;
    assert_eq!(
      10, effective,
      "effective length should exclude union of parent+child non_char positions"
    );

    // Non-char positions (N, gap) are tracked through non_char ranges, not
    // Fitch substitutions. Applying child-side subs+indels to root composition
    // does NOT reproduce child composition when the child has N or gaps that
    // the root does not.
    //
    // Child A has N at positions 8-9 where root has G,T. No Fitch sub exists
    // for these positions because N is ambiguous (non-informative).
    let child_node = &sparse.nodes[&a_key];
    let child_edge = &sparse.edges[&child_side_key];

    let mut derived = root_node.seq.composition.clone();
    for sub in child_edge.fitch_subs() {
      derived.add_sub(sub);
    }
    for indel in &child_edge.indels {
      derived.add_indel(indel);
    }
    assert_ne!(
      derived, child_node.seq.composition,
      "root + subs/indels should NOT equal child composition (non-char positions are not Fitch subs)"
    );

    Ok(())
  }

  /// Reroot on branch AB->A with trivial root removal (edge_merge path).
  ///
  /// Same tree and alignment as `test_fitch_reroot_sparse_on_branch_ab_to_a`.
  /// After splitting AB->A and inverting the path to the old root, the old root
  /// becomes trivial (1 parent, 1 child). `remove_node_if_trivial` merges the
  /// two edges around it. Verifies:
  ///
  /// - Root sequence is derived from inverted edges (not from merge parent_edge)
  /// - Old root node is removed from the partition
  /// - Merged edge has composed substitutions from both original edges
  /// - New root node fields (gaps, non_char) are correct
  #[test]
  fn test_fitch_reroot_sparse_with_trivial_root_removal() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let alphabet = Alphabet::default();

    let partitions = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions, &aln)?;

    let fitch = Arc::try_unwrap(partitions.into_iter().next().unwrap())
      .map(|rw| rw.into_inner())
      .map_err(|_e| eyre::eyre!("Fitch partition still shared"))?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let mut sparse = fitch.into_marginal_sparse(gtr, &graph)?;

    let old_root_key = graph.get_exactly_one_root()?.read_arc().key();
    let ab_key = find_node_key_by_name(&graph, "AB").expect("AB node not found");
    let a_key = find_node_key_by_name(&graph, "A").expect("A node not found");

    let edge_ab_a_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == ab_key && e.target() == a_key
      })
      .map(|e| e.read_arc().key())
      .expect("AB->A edge not found");

    // Record original root->CD subs (will appear in merged edge)
    let edge_root_cd_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        let cd_key = find_node_key_by_name(&graph, "CD").unwrap();
        e.source() == old_root_key && e.target() == cd_key
      })
      .map(|e| e.read_arc().key())
      .expect("root->CD edge not found");

    let orig_root_cd_subs: Vec<String> = sparse.edges[&edge_root_cd_key]
      .fitch_subs()
      .iter()
      .map(|s| s.to_string())
      .collect();

    let orig_root_ab_subs: Vec<String> = {
      let edge_root_ab_key = graph
        .get_edges()
        .iter()
        .find(|e| {
          let e = e.read_arc();
          e.source() == old_root_key && e.target() == ab_key
        })
        .map(|e| e.read_arc().key())
        .expect("root->AB edge not found");
      sparse.edges[&edge_root_ab_key]
        .fitch_subs()
        .iter()
        .map(|s| s.to_string())
        .collect()
    };

    // Split AB->A at midpoint
    let split_info = split_edge(&mut graph, edge_ab_a_key, 0.5)?;
    let new_root_key = split_info.new_node_key;

    // Invert path from new_root to old_root
    let inverted_edge_keys = apply_reroot_topology(&mut graph, old_root_key, new_root_key)?;

    // Old root is now trivial (1 parent: AB, 1 child: CD) - remove it
    let edge_merge = remove_node_if_trivial(&mut graph, old_root_key)?;
    assert!(edge_merge.is_some(), "old root should be trivial after reroot");

    let changes = RerootChanges {
      edge_split: Some(split_info),
      edge_merge,
      inverted_edge_keys,
    };

    sparse.apply_reroot(&changes)?;

    // Root sequence should be AB's ancestral state (derived from inverted edges)
    assert_eq!(
      "ACATCCCTGTA--G--",
      sparse.root_sequence.as_str(),
      "root_sequence should equal AB's sequence after reroot with merge"
    );

    // Old root node should be removed
    assert!(
      !sparse.nodes.contains_key(&old_root_key),
      "old root node should be removed from partition after trivial root removal"
    );

    // New root node should have correct fields
    let root_node = &sparse.nodes[&new_root_key];
    assert_eq!(vec![(11, 13), (14, 16)], root_node.seq.gaps);
    assert!(root_node.seq.unknown.is_empty());
    assert_eq!(root_node.seq.gaps, root_node.seq.non_char);

    // Merged edge (AB->CD, replacing AB->root + root->CD) should have composed subs.
    // Original root->AB: [G4T, A7C], inverted to AB->root: [T4G, C7A]
    // Original root->CD: [A1C, A3G]
    // After merge removal, the merged edge connects AB->CD.
    // Composed subs: AB->root subs [T4G, C7A] chained with root->CD subs [A1C, A3G]
    let merge_info = changes.edge_merge.as_ref().unwrap();
    let merged_edge = &sparse.edges[&merge_info.merged_edge_key];
    assert!(
      !merged_edge.fitch_subs().is_empty(),
      "merged edge should have composed substitutions from both original edges"
    );

    // Sanity: original subs are what we expect
    assert_eq!(vec_of_owned!["G4T", "A7C"], orig_root_ab_subs);
    assert_eq!(vec_of_owned!["A1C", "A3G"], orig_root_cd_subs);

    Ok(())
  }

  /// After reroot, a marginal forward pass should produce non-zero fixed_counts.
  ///
  /// This is the integration test that would have caught the original bug where
  /// `SparseNodePartition::empty()` left composition as zero, causing
  /// `process_node_forward` to seed `msg_to_child.fixed_counts` with zero
  /// multiplicity for all characters.
  #[test]
  fn test_fitch_reroot_sparse_forward_pass_nonzero_fixed_counts() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
      &*NUC_ALPHABET,
    )?;

    let mut graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let alphabet = Alphabet::default();

    let partitions_fitch = [Arc::new(RwLock::new(PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions_fitch, &aln)?;

    let fitch = Arc::try_unwrap(partitions_fitch.into_iter().next().unwrap())
      .map(|rw| rw.into_inner())
      .map_err(|_e| eyre::eyre!("Fitch partition still shared"))?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let sparse = fitch.into_marginal_sparse(gtr, &graph)?;

    // Run initial marginal pass before reroot
    let partitions: Vec<Arc<RwLock<PartitionMarginalSparse>>> = vec![Arc::new(RwLock::new(sparse))];
    update_marginal(&graph, &partitions)?;

    // Reroot on AB->A
    let old_root_key = graph.get_exactly_one_root()?.read_arc().key();
    let ab_key = find_node_key_by_name(&graph, "AB").expect("AB node not found");
    let a_key = find_node_key_by_name(&graph, "A").expect("A node not found");

    let edge_ab_a_key = graph
      .get_edges()
      .iter()
      .find(|e| {
        let e = e.read_arc();
        e.source() == ab_key && e.target() == a_key
      })
      .map(|e| e.read_arc().key())
      .expect("AB->A edge not found");

    let split_info = split_edge(&mut graph, edge_ab_a_key, 0.5)?;
    let new_root_key = split_info.new_node_key;
    let inverted_edge_keys = apply_reroot_topology(&mut graph, old_root_key, new_root_key)?;

    let changes = RerootChanges {
      edge_split: Some(split_info),
      edge_merge: None,
      inverted_edge_keys,
    };

    partitions[0].write_arc().apply_reroot(&changes)?;

    // Run marginal pass after reroot
    update_marginal(&graph, &partitions)?;

    // Check that msg_to_child.fixed_counts for edges from the new root have
    // non-zero entries. With zero composition (the original bug), all counts
    // would be zero.
    let partition = partitions[0].read_arc();
    for edge_ref in graph.get_edges() {
      let edge = edge_ref.read_arc();
      if edge.source() == new_root_key {
        let edge_data = &partition.edges[&edge.key()];
        let total: usize = edge_data.msg_to_child.fixed_counts.counts().values().sum();
        assert!(
          total > 0,
          "msg_to_child.fixed_counts for edge from new root to {:?} should have non-zero entries (total={total})",
          edge.target()
        );
      }
    }

    Ok(())
  }
}
