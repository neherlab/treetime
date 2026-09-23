#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::compress_sequences;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};

  use crate::partition::fitch::partition::PartitionFitch;
  use crate::partition::marginal::sparse::reroot::reroot_sparse;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::composition::Composition;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  use treetime_graph::graph::Graph;

  use treetime_graph::reroot::{
    RerootChanges, apply_reroot_topology, record_split, remove_node_if_trivial, split_edge, trivial_node_branch_lengths,
  };
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AlignmentRecord, AsciiChar};

  use crate::ancestral::__tests__::test_fitch::tests::helpers::*;
  use treetime_utils::vec_of_owned;

  #[test]
  fn test_fitch_reroot_sparse_on_branch_ab_to_a() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let alphabet = Alphabet::default();

    let mut fitch = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut fitch, &node_seq_inputs(&graph, &names, aln))?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, gtr, node_states);

    let old_root_key = graph.get_exactly_one_root()?.key();
    let ab_key = find_node_key_by_name(&graph, &names, "AB").expect("AB node not found");
    let a_key = find_node_key_by_name(&graph, &names, "A").expect("A node not found");

    let edge_ab_a_key = graph
      .get_edges()
      .find(|e| e.source() == ab_key && e.target() == a_key)
      .map(|e| e.key())
      .expect("AB->A edge not found");

    let orig_subs: Vec<String> = recon.partition.obs_edges[&edge_ab_a_key]
      .fitch_subs()
      .iter()
      .map(|s| s.to_string())
      .collect();
    let orig_indels: Vec<String> = recon.partition.obs_edges[&edge_ab_a_key]
      .indels
      .iter()
      .map(|i| i.to_string())
      .collect();

    let edge_root_ab_key = graph
      .get_edges()
      .find(|e| e.source() == old_root_key && e.target() == ab_key)
      .map(|e| e.key())
      .expect("root->AB edge not found");

    let orig_root_ab_subs: Vec<String> = recon.partition.obs_edges[&edge_root_ab_key]
      .fitch_subs()
      .iter()
      .map(|s| s.to_string())
      .collect();
    let orig_root_ab_indels: Vec<String> = recon.partition.obs_edges[&edge_root_ab_key]
      .indels
      .iter()
      .map(|i| i.to_string())
      .collect();

    let split_info = split_edge(&mut graph, edge_ab_a_key, 0.5, branch_lengths[&edge_ab_a_key])?;
    let new_root_key = split_info.new_node_key;
    let parent_side_key = split_info.parent_side_edge_key;
    let child_side_key = split_info.child_side_edge_key;

    let inverted_edge_keys = apply_reroot_topology(&mut graph, old_root_key, new_root_key)?;

    let changes = RerootChanges {
      edge_split: Some(split_info),
      edge_merge: None,
      inverted_edge_keys,
    };

    let recon = reroot_sparse(recon.partition, recon.gtr, recon.node_states, &changes)?;

    let expected_root_seq = "ACATCCCTGTA--G--";
    let actual_root_seq = recon.partition.root_sequence.as_str();
    assert_eq!(
      expected_root_seq, actual_root_seq,
      "root_sequence after reroot should equal AB's sequence"
    );

    let child_edge = &recon.partition.obs_edges[&child_side_key];
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

    let parent_edge = &recon.partition.obs_edges[&parent_side_key];
    assert!(
      parent_edge.fitch_subs().is_empty(),
      "parent-side edge should have no subs"
    );
    assert!(parent_edge.indels.is_empty(), "parent-side edge should have no indels");

    let inv_edge = &recon.partition.obs_edges[&edge_root_ab_key];
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

    assert_eq!(vec_of_owned!["G4T", "A7C"], orig_root_ab_subs);
    assert_eq!(vec_of_owned!["11--13: TT -> --"], orig_root_ab_indels);

    let root_node = &recon.partition.obs_nodes[&new_root_key];
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
      expected_comp, root_node.composition,
      "new root node composition should match root_sequence character counts"
    );

    assert_eq!(
      vec![(11, 13), (14, 16)],
      root_node.gaps,
      "new root gaps should cover gap positions in root_sequence"
    );
    assert!(
      root_node.unknown.is_empty(),
      "new root should have no unknown (N) positions"
    );
    assert_eq!(
      root_node.gaps, root_node.non_char,
      "non_char should equal gaps when there are no N positions"
    );

    let effective = recon.partition.edge_effective_length(&graph, child_side_key)?;
    assert_eq!(
      10, effective,
      "effective length should exclude union of parent+child non_char positions"
    );

    let child_node = &recon.partition.obs_nodes[&a_key];
    let child_edge = &recon.partition.obs_edges[&child_side_key];

    let mut derived = root_node.composition.clone();
    for sub in child_edge.fitch_subs() {
      derived.add_sub(sub);
    }
    for indel in &child_edge.indels {
      derived.add_indel(indel);
    }
    assert_ne!(
      derived, child_node.composition,
      "root + subs/indels should NOT equal child composition (non-char positions are not Fitch subs)"
    );

    Ok(())
  }

  #[test]
  fn test_fitch_reroot_sparse_with_trivial_root_removal() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let alphabet = Alphabet::default();

    let mut fitch = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut fitch, &node_seq_inputs(&graph, &names, aln))?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, gtr, node_states);

    let old_root_key = graph.get_exactly_one_root()?.key();
    let ab_key = find_node_key_by_name(&graph, &names, "AB").expect("AB node not found");
    let a_key = find_node_key_by_name(&graph, &names, "A").expect("A node not found");

    let edge_ab_a_key = graph
      .get_edges()
      .find(|e| e.source() == ab_key && e.target() == a_key)
      .map(|e| e.key())
      .expect("AB->A edge not found");

    let edge_root_cd_key = graph
      .get_edges()
      .find(|e| {
        let cd_key = find_node_key_by_name(&graph, &names, "CD").unwrap();
        e.source() == old_root_key && e.target() == cd_key
      })
      .map(|e| e.key())
      .expect("root->CD edge not found");

    let orig_root_cd_subs: Vec<String> = recon.partition.obs_edges[&edge_root_cd_key]
      .fitch_subs()
      .iter()
      .map(|s| s.to_string())
      .collect();

    let orig_root_ab_subs: Vec<String> = {
      let edge_root_ab_key = graph
        .get_edges()
        .find(|e| e.source() == old_root_key && e.target() == ab_key)
        .map(|e| e.key())
        .expect("root->AB edge not found");
      recon.partition.obs_edges[&edge_root_ab_key]
        .fitch_subs()
        .iter()
        .map(|s| s.to_string())
        .collect()
    };

    let split_info = split_edge(&mut graph, edge_ab_a_key, 0.5, branch_lengths[&edge_ab_a_key])?;
    let new_root_key = split_info.new_node_key;

    let inverted_edge_keys = apply_reroot_topology(&mut graph, old_root_key, new_root_key)?;

    let (old_root_parent, old_root_child) = trivial_node_branch_lengths(&graph, old_root_key, &branch_lengths);
    let edge_merge = remove_node_if_trivial(&mut graph, old_root_key, old_root_parent, old_root_child)?;
    assert!(edge_merge.is_some(), "old root should be trivial after reroot");

    let changes = RerootChanges {
      edge_split: Some(split_info),
      edge_merge,
      inverted_edge_keys,
    };

    let recon = reroot_sparse(recon.partition, recon.gtr, recon.node_states, &changes)?;

    assert_eq!(
      "ACATCCCTGTA--G--",
      recon.partition.root_sequence.as_str(),
      "root_sequence should equal AB's sequence after reroot with merge"
    );

    assert!(
      !recon.partition.obs_nodes.contains_key(&old_root_key),
      "old root node should be removed from partition after trivial root removal"
    );

    let root_node = &recon.partition.obs_nodes[&new_root_key];
    assert_eq!(vec![(11, 13), (14, 16)], root_node.gaps);
    assert!(root_node.unknown.is_empty());
    assert_eq!(root_node.gaps, root_node.non_char);

    let merge_info = changes.edge_merge.as_ref().unwrap();
    let merged_edge = &recon.partition.obs_edges[&merge_info.merged_edge_key];

    assert_eq!(vec_of_owned!["G4T", "A7C"], orig_root_ab_subs);
    assert_eq!(vec_of_owned!["A1C", "A3G"], orig_root_cd_subs);

    let mut merged_subs: Vec<String> = merged_edge.fitch_subs().iter().map(|s| s.to_string()).collect();
    merged_subs.sort();
    assert_eq!(vec_of_owned!["A1C", "A3G", "C7A", "T4G"], merged_subs);

    Ok(())
  }

  #[test]
  fn test_fitch_reroot_sparse_forward_pass_nonzero_fixed_counts() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let alphabet = Alphabet::default();

    let mut fitch = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut fitch, &node_seq_inputs(&graph, &names, aln))?;

    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      ..JC69Params::default()
    })?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, gtr, node_states);

    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let old_root_key = graph.get_exactly_one_root()?.key();
    let ab_key = find_node_key_by_name(&graph, &names, "AB").expect("AB node not found");
    let a_key = find_node_key_by_name(&graph, &names, "A").expect("A node not found");

    let edge_ab_a_key = graph
      .get_edges()
      .find(|e| e.source() == ab_key && e.target() == a_key)
      .map(|e| e.key())
      .expect("AB->A edge not found");

    let mut branch_lengths = branch_lengths;
    let split_info = split_edge(&mut graph, edge_ab_a_key, 0.5, branch_lengths[&edge_ab_a_key])?;
    record_split(&mut branch_lengths, &split_info);
    let new_root_key = split_info.new_node_key;
    let inverted_edge_keys = apply_reroot_topology(&mut graph, old_root_key, new_root_key)?;

    let changes = RerootChanges {
      edge_split: Some(split_info),
      edge_merge: None,
      inverted_edge_keys,
    };

    let recon = reroot_sparse(recon.partition, recon.gtr, recon.node_states, &changes)?;

    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let root_edge_totals: Vec<(_, usize)> = graph
      .get_edges()
      .filter_map(|edge_ref| {
        let edge = edge_ref;
        (edge.source() == new_root_key).then(|| {
          let edge_data = &recon.edges.forward[&edge.key()];
          let total: usize = edge_data.msg_to_child.fixed_counts.counts().values().sum();
          (edge.target(), total)
        })
      })
      .collect();
    assert!(!root_edge_totals.is_empty(), "new root should have outgoing edges");
    for (target, total) in &root_edge_totals {
      assert!(
        *total > 0,
        "msg_to_child.fixed_counts for edge from new root to {target:?} should have non-zero entries (total={total})"
      );
    }

    Ok(())
  }
}
