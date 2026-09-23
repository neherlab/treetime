#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::fitch::{
    ancestral_reconstruction_fitch, attach_seqs_to_graph, compress_sequences, fitch_backward, fitch_forward,
  };

  use crate::o;
  use crate::partition::fitch::partition::PartitionFitch;

  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;

  use eyre::Report;
  use indoc::indoc;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::LazyLock;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;

  use helpers::*;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::io::json::{JsonPretty, json_write_str};
  use treetime_utils::vec_of_owned;

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

  pub mod helpers {
    use super::*;

    pub fn get_node_name(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> String {
      names[&key].clone().expect("node has name")
    }

    pub fn collect_edge_subs(
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

    pub fn get_root_variable_positions(graph: &Graph, partition: &PartitionFitch) -> Vec<usize> {
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

    pub fn get_root_state_sets(graph: &Graph, partition: &PartitionFitch) -> BTreeMap<usize, String> {
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

    pub fn get_node_variable_positions_by_name(
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

    pub fn get_root_seq(graph: &Graph, partition: &PartitionFitch) -> String {
      let root = graph.get_exactly_one_root().expect("graph has exactly one root");
      let root_key = root.key();
      partition.nodes[&root_key].seq.sequence.as_str().to_owned()
    }

    pub fn get_internal_sequences(
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

    pub fn collect_edge_indels(
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

    pub static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);
  }
}
