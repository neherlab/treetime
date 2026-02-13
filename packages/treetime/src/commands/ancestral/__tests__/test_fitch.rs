#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{
    ancestral_reconstruction_fitch, attach_seqs_to_graph, compress_sequences, fitch_backward, fitch_forward,
    get_common_length,
  };
  use crate::o;
  use crate::representation::graph_ancestral::GraphAncestral;
  use crate::representation::partition_fitch::PartitionFitch;
  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::json::{JsonPretty, json_write_str};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::vec_of_owned;

  fn get_node_name(graph: &GraphAncestral, key: treetime_graph::node::GraphNodeKey) -> String {
    let node = graph.get_node(key).expect("node exists");
    node
      .read_arc()
      .payload()
      .read_arc()
      .name
      .clone()
      .expect("node has name")
  }

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
          .subs
          .iter()
          .map(|sub| sub.to_string())
          .collect_vec();
        (edge_name, subs)
      })
      .collect()
  }

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

  fn get_root_seq(graph: &GraphAncestral, partition: &PartitionFitch) -> String {
    let root = graph.get_exactly_one_root().expect("graph has exactly one root");
    let root_key = root.read_arc().key();
    partition.nodes[&root_key].seq.sequence.as_str().to_owned()
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

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  #[test]
  fn test_ancestral_reconstruction_fitch() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

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
        ACATCCCTGTA-TG--
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
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_ancestral_reconstruction_fitch_with_leaves() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

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
        ACATCCCTGTA-TG--
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
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_fitch_internals() -> Result<(), Report> {
    drop(rayon::ThreadPoolBuilder::new().num_threads(1).build_global());

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
      o!("AB->A")     => vec_of_owned!["12--13: T -> -", "14--16: -- -> AC"],
      o!("AB->B")     => vec![],
      o!("root->AB")  => vec_of_owned!["11--12: T -> -"],
      o!("CD->C")     => vec![],
      o!("CD->D")     => vec![],
      o!("root->CD")  => vec![],
    };
    assert_eq!(expected_indels, actual_indels);

    Ok(())
  }

  #[test]
  fn test_fitch_complex_gaps() -> Result<(), Report> {
    drop(rayon::ThreadPoolBuilder::new().num_threads(1).build_global());

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

  #[test]
  fn test_fitch_polytomy() -> Result<(), Report> {
    drop(rayon::ThreadPoolBuilder::new().num_threads(1).build_global());

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
      o!("CDE->D")    => vec![o!("C3T"), o!("C4T")],
      o!("CDE->E")    => vec![],
      o!("root->CDE") => vec![o!("A4C"), o!("C2G")],
    };
    assert_eq!(expected_subs, actual_subs);

    // Verify indels on edges
    let actual_indels = collect_edge_indels(&graph, &partition);
    let expected_indels = btreemap! {
      o!("AB->A")     => vec![o!("3--4: A -> -")],
      o!("AB->B")     => vec![o!("1--2: C -> -")],
      o!("root->AB")  => vec![],
      o!("CDE->C")    => vec![o!("2--3: C -> -")],
      o!("CDE->D")    => vec![],
      o!("CDE->E")    => vec![],
      o!("root->CDE") => vec![o!("2--3: - -> C")],
    };
    assert_eq!(expected_indels, actual_indels);

    Ok(())
  }

  #[test]
  fn test_fitch_backward_state() -> Result<(), Report> {
    drop(rayon::ThreadPoolBuilder::new().num_threads(1).build_global());

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
          0 => o!("{A, C, G, T}"),  // union of all: A(leaf A), G(leaf B), C(leaf C), T(leaf D)
          2 => o!("{A, G}"),        // A from AB subtree; G from CD subtree
          3 => o!("{G, T}"),        // T from AB subtree; G from CD subtree
          5 => o!("{C, G}"),        // C from AB subtree; G from CD subtree
          6 => o!("{A, C, G}"),     // union: C from AB, then A or G from CD
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
    fitch_forward(&graph, &partitions);

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
}
