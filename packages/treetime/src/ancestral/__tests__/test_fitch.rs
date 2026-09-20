#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::{
    ancestral_reconstruction_fitch, attach_seqs_to_graph, compress_sequences, fitch_backward, fitch_forward,
  };
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::o;
  use crate::partition::fitch::partition::PartitionFitch;
  use crate::partition::marginal::sparse::reroot::reroot_sparse;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::composition::Composition;
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::LazyLock;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_graph::reroot::{
    RerootChanges, apply_reroot_topology, record_split, remove_node_if_trivial, split_edge, trivial_node_branch_lengths,
  };
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AlignmentRecord, AsciiChar};
  use treetime_utils::io::json::{JsonPretty, json_write_str};
  use treetime_utils::vec_of_owned;

  fn get_node_name(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> String {
    names[&key].clone().expect("node has name")
  }

  fn collect_edge_subs(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    partition: &PartitionFitch,
  ) -> BTreeMap<String, Vec<String>> {
    graph
      .get_edges()
      .map(|edge| {
        let parent_name = get_node_name(names, edge.source());
        let child_name = get_node_name(names, edge.target());
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

  fn get_root_variable_positions(graph: &Graph, partition: &PartitionFitch) -> Vec<usize> {
    let root = graph.get_exactly_one_root().expect("graph has exactly one root");
    let root_key = root.key();
    partition.nodes[&root_key]
      .seq
      .fitch
      .variable
      .keys()
      .copied()
      .collect_vec()
  }

  fn get_root_state_sets(graph: &Graph, partition: &PartitionFitch) -> BTreeMap<usize, String> {
    let root = graph.get_exactly_one_root().expect("graph has exactly one root");
    let root_key = root.key();
    partition.nodes[&root_key]
      .seq
      .fitch
      .variable
      .iter()
      .map(|(&pos, states)| (pos, states.to_string()))
      .collect()
  }

  fn get_node_variable_positions_by_name(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    partition: &PartitionFitch,
    name: &str,
  ) -> Vec<usize> {
    for node in graph.get_nodes() {
      let node_name = names[&node.key()].clone();
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

  fn get_root_seq(graph: &Graph, partition: &PartitionFitch) -> String {
    let root = graph.get_exactly_one_root().expect("graph has exactly one root");
    let root_key = root.key();
    partition.nodes[&root_key].seq.sequence.as_str().to_owned()
  }

  fn get_internal_sequences(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    partition: &PartitionFitch,
  ) -> BTreeMap<String, String> {
    graph
      .get_internal_nodes()
      .map(|node| {
        let node_name = names[&node.key()].clone().unwrap();
        let sequence = partition.nodes[&node.key()].seq.sequence.as_str().to_owned();
        (node_name, sequence)
      })
      .collect()
  }

  fn collect_edge_indels(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    partition: &PartitionFitch,
  ) -> BTreeMap<String, Vec<String>> {
    graph
      .get_edges()
      .map(|edge| {
        let parent_name = get_node_name(names, edge.source());
        let child_name = get_node_name(names, edge.target());
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

  #[test]
  fn test_ancestral_reconstruction_fitch() -> Result<(), Report> {
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

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let alphabet = Alphabet::default();

    let mut partition = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut partition, &node_seq_inputs(&graph, &names, aln))?;
    let mut partitions_parsimony = [partition];

    let mut actual = BTreeMap::new();
    let emitted = ancestral_reconstruction_fitch(&graph, false, &mut partitions_parsimony)?;
    for key in emitted {
      actual.insert(
        names[&key].clone(),
        partitions_parsimony[0].node_sequence(key).to_string(),
      );
    }

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_ancestral_reconstruction_fitch_with_leaves() -> Result<(), Report> {
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

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let alphabet = Alphabet::default();

    let mut partition = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut partition, &node_seq_inputs(&graph, &names, aln))?;
    let mut partitions_parsimony = [partition];

    let mut actual = BTreeMap::new();
    let emitted = ancestral_reconstruction_fitch(&graph, true, &mut partitions_parsimony)?;
    for key in emitted {
      actual.insert(
        names[&key].clone(),
        partitions_parsimony[0].node_sequence(key).to_string(),
      );
    }

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_compress_sequences_retains_internal_exact_sequences() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let alphabet = Alphabet::default();

    let mut partition = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut partition, &node_seq_inputs(&graph, &names, aln))?;
    let actual = get_internal_sequences(&graph, &names, &partition);
    let expected = btreemap! {
      o!("AB") => o!("GCGTACGT"),
      o!("CD") => o!("GCGTACGT"),
      o!("root") => o!("GCGTACGT"),
    };

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_fitch_internals() -> Result<(), Report> {
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
    let graph: Graph = graph;

    let alphabet = Alphabet::default();

    let mut partition = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut partition, &node_seq_inputs(&graph, &names, aln))?;

    let actual_subs = collect_edge_subs(&graph, &names, &partition);
    let expected_subs = btreemap! {
      o!("AB->A")    => vec_of_owned!["C6G", "T8C"],
      o!("AB->B")    => vec_of_owned!["A1G"],
      o!("root->AB") => vec_of_owned!["G4T", "A7C"],
      o!("CD->C")    => vec_of_owned!["C6G"],
      o!("CD->D")    => vec_of_owned!["C1T", "A7G"],
      o!("root->CD") => vec_of_owned!["A1C", "A3G"],
    };
    assert_eq!(expected_subs, actual_subs);

    let actual_indels = collect_edge_indels(&graph, &names, &partition);
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

  #[test]
  fn test_fitch_complex_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,(D:0.05,E:0.03)DE:0.01)CDE:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let alphabet = Alphabet::default();

    let mut partition = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut partition, &node_seq_inputs(&graph, &names, aln))?;

    let actual_subs = collect_edge_subs(&graph, &names, &partition);
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

    let actual_indels = collect_edge_indels(&graph, &names, &partition);
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

  #[test]
  fn test_fitch_polytomy() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.05,E:0.03)CDE:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let alphabet = Alphabet::default();

    let mut partition = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };
    compress_sequences(&graph, &mut partition, &node_seq_inputs(&graph, &names, aln))?;

    let actual_subs = collect_edge_subs(&graph, &names, &partition);
    let expected_subs = btreemap! {
      o!("AB->A")     => vec![],
      o!("AB->B")     => vec![],
      o!("root->AB")  => vec![],
      o!("CDE->C")    => vec![],
      o!("CDE->D")    => vec![],
      o!("CDE->E")    => vec![o!("T4C")],
      o!("root->CDE") => vec![o!("C2G"), o!("A4T")],
    };
    assert_eq!(expected_subs, actual_subs);

    let actual_indels = collect_edge_indels(&graph, &names, &partition);
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

  #[test]
  fn test_fitch_backward_state() -> Result<(), Report> {
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
    let graph: Graph = graph;

    let alphabet = Alphabet::default();

    let mut partition = PartitionFitch {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    attach_seqs_to_graph(&graph, &mut partition, &node_seq_inputs(&graph, &names, aln))?;
    fitch_backward(&graph, &mut partition)?;

    {
      let variable_positions = get_root_variable_positions(&graph, &partition);
      assert_eq!(vec![0, 2, 3, 5, 6], variable_positions);

      let root_seq = get_root_seq(&graph, &partition);
      assert_eq!("~C~~C~~TGTATTG..", root_seq);

      let state_sets = get_root_state_sets(&graph, &partition);
      assert_eq!(
        btreemap! {
          0 => o!("{A, C, G, T}"),
          2 => o!("{A, G}"),
          3 => o!("{G, T}"),
          5 => o!("{C, G}"),
          6 => o!("{A, C, G}"),
        },
        state_sets
      );

      let ab_vars = get_node_variable_positions_by_name(&graph, &names, &partition, "AB");
      assert!(!ab_vars.is_empty(), "AB should have variable positions");

      let cd_vars = get_node_variable_positions_by_name(&graph, &names, &partition, "CD");
      assert!(!cd_vars.is_empty(), "CD should have variable positions");
    }

    fitch_forward(&graph, &mut partition)?;

    {
      let variable_positions = get_root_variable_positions(&graph, &partition);
      assert_eq!(vec![0, 2, 3, 5, 6], variable_positions);

      let root_seq = get_root_seq(&graph, &partition);
      assert_eq!("ACAGCCATGTATTG--", root_seq);
    }

    Ok(())
  }

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
