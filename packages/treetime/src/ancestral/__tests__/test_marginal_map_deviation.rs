#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
  use crate::partition::marginal::sample::SampleMode;
  use crate::seq::alignment::node_seq_inputs;
  use crate::test_utils::{emitted_sequences_by_name, node_keys};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::sync::random::get_random_number_generator;

  #[test]
  fn test_marginal_map_deviation_does_not_leak_into_subtree() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >T1
      TAAAAAAAAA
      >C1
      CAAAAAAAAA
      >C2
      CAAAAAAAAA
      >C3
      CAAAAAAAAA
      >C4
      CAAAAAAAAA
      >C5
      CAAAAAAAAA
      >C6
      CAAAAAAAAA
      >C7
      CAAAAAAAAA
    "#})?;
    let nwk_parsed = nwk_read(
      b"((T1:0.0,((((C1:0.005,C2:0.005)Y3:0.005,C3:0.005)Y2:0.005,C4:0.005)Y1:0.005,C5:0.005)Z:0.5,C6:0.5)X:0.3,C7:0.3)root:0.0;".as_slice(),
    )?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let sparse = reconstruct_sparse(&graph, &branch_lengths, &names, &aln)?;
    let dense = reconstruct_dense(&graph, &branch_lengths, &names, &aln)?;

    assert_eq!('T', nuc_at(&sparse, "X", 0), "X is the node whose argmax deviates");

    for node in ["Z", "Y1", "Y2", "Y3"] {
      assert_eq!(
        'C',
        nuc_at(&sparse, node, 0),
        "{node} resolved position 0 to C and must not inherit X's T"
      );
    }

    assert_eq!(dense, sparse, "sparse must reproduce the dense reconstruction");
    Ok(())
  }

  #[test]
  fn test_marginal_map_deviation_keeps_inherited_deletions() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >D1
      -CGTACGTAC
      >D2
      -CGTACGTAC
      >D3
      -CGTACGTAC
      >A1
      ACGTACGTAC
      >A2
      TCGTACGTAC
    "#})?;
    let nwk_parsed =
      nwk_read(b"(((D1:0.05,D2:0.05)DD:0.05,D3:0.05)DEL:0.2,(A1:0.05,A2:0.05)POLY:0.2)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let sparse = reconstruct_sparse(&graph, &branch_lengths, &names, &aln)?;

    for node in ["DEL", "DD"] {
      assert_eq!(
        '-',
        nuc_at(&sparse, node, 0),
        "{node} is deleted at position 0 and must not have a residue restored"
      );
    }

    assert_eq!(
      reconstruct_dense(&graph, &branch_lengths, &names, &aln)?,
      sparse,
      "sparse must reproduce the dense reconstruction"
    );
    Ok(())
  }

  fn nuc_at(seqs: &BTreeMap<String, String>, node: &str, pos: usize) -> char {
    seqs[node].chars().nth(pos).expect("position within sequence")
  }

  fn parse_aln(fasta: &str) -> Result<Vec<AlignmentRecord>, Report> {
    Ok(
      fasta_read(fasta.as_bytes(), &Alphabet::default())?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect(),
    )
  }

  fn reconstruct_sparse(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
  ) -> Result<BTreeMap<String, String>, Report> {
    let fitch = create_fitch_partition(
      graph,
      0,
      Alphabet::default(),
      node_seq_inputs(graph, names, aln.to_vec()),
    )?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let recon = MarginalReconstruction::Sparse(SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    ));
    let (recon, _) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    named_strings(
      names,
      &node_keys(graph),
      &recon.sample_sequences(graph, SampleMode::Argmax, &mut get_random_number_generator(0))?,
      |key| recon.node_sequence(graph, false, key),
    )
  }

  fn reconstruct_dense(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
  ) -> Result<BTreeMap<String, String>, Report> {
    let partition = PartitionMarginalDense::new(
      0,
      Alphabet::default(),
      graph,
      &node_seq_inputs(graph, names, aln.to_vec()),
    )?;
    let recon = MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, jc69(JC69Params::default())?));
    let (recon, _) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    named_strings(
      names,
      &node_keys(graph),
      &recon.sample_sequences(graph, SampleMode::Argmax, &mut get_random_number_generator(0))?,
      |key| recon.node_sequence(graph, false, key),
    )
  }

  fn named_strings(
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    emitted_nodes: &[GraphNodeKey],
    sampled: &BTreeMap<GraphNodeKey, Seq>,
    node_sequence: impl Fn(GraphNodeKey) -> Result<Seq, Report>,
  ) -> Result<BTreeMap<String, String>, Report> {
    Ok(
      emitted_sequences_by_name(names, emitted_nodes, sampled, node_sequence)?
        .into_iter()
        .map(|(name, seq)| (name, seq.as_str().to_owned()))
        .collect(),
    )
  }
}
