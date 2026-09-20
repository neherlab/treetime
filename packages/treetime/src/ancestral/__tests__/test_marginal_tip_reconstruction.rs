#![allow(
  clippy::disallowed_methods,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::{ancestral_reconstruction, branch_lengths_or_zero};
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::ancestral::sample::SampleMode;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AlignmentRecord, Seq};

  #[test]
  fn test_marginal_tip_c3_preserves_observed_leaf() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >A
      ACGT
      >B
      GCGT
    "#})?;
    let nwk_parsed = nwk_read_str("(A:0.4,B:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let sparse = reconstruct_sparse(&graph, &branch_lengths, &names, &aln, false)?;
    let dense = reconstruct_dense(&graph, &branch_lengths, &names, &aln, false)?;

    assert_eq!("ACGT", sparse["A"]);
    assert_eq!("GCGT", sparse["B"]);
    assert_eq!(sparse, dense);
    Ok(())
  }

  #[test]
  fn test_marginal_tip_parent_plus_muts_equals_child() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >A
      ACGTA
      >B
      AGGTA
      >C
      ATGTC
    "#})?;
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)AB:0.1,C:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let alphabet = Alphabet::default();
    let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = SparseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states);
    let (mut recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let seqs = reconstruct_named_sparse(&graph, &names, &mut recon, false)?;

    for edge in graph.get_edges() {
      let parent = node_name(&names, edge.source());
      let child = node_name(&names, edge.target());
      let mut expected = seqs[&parent].clone();
      for sub in recon.edge_subs(edge.key())? {
        expected[sub.pos()] = sub.qry();
      }
      assert_eq!(
        seqs[&child].as_str().to_owned(),
        expected.as_str().to_owned(),
        "parent {parent} + muts must equal child {child}"
      );
    }
    Ok(())
  }

  #[test]
  fn test_marginal_tip_impute_resolves_n_and_iupac() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >A
      ANRT
      >B
      ACGT
      >C
      ACGT
    "#})?;
    let nwk_parsed = nwk_read_str("(A:0.1,(B:0.1,C:0.1)BC:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let sparse_plain = reconstruct_sparse(&graph, &branch_lengths, &names, &aln, false)?;
    let dense_plain = reconstruct_dense(&graph, &branch_lengths, &names, &aln, false)?;
    assert_eq!("ANRT", sparse_plain["A"]);
    assert_eq!(sparse_plain, dense_plain);

    let sparse_imputed = reconstruct_sparse(&graph, &branch_lengths, &names, &aln, true)?;
    let dense_imputed = reconstruct_dense(&graph, &branch_lengths, &names, &aln, true)?;
    assert_eq!("ACGT", sparse_imputed["A"]);
    assert_eq!(sparse_imputed, dense_imputed);
    Ok(())
  }

  #[test]
  fn test_gm_marginal_tip_impute_matches_v0() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >A
      ANRT
      >B
      ACGT
      >C
      ACGT
    "#})?;
    let nwk_parsed = nwk_read_str("(A:0.1,(B:0.1,C:0.1)BC:0.1)root:0.0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let sparse = reconstruct_sparse(&graph, &branch_lengths, &names, &aln, true)?;
    let dense = reconstruct_dense(&graph, &branch_lengths, &names, &aln, true)?;

    assert_eq!("ACGT", sparse["A"]);
    assert_eq!("ACGT", dense["A"]);
    Ok(())
  }

  fn parse_aln(fasta: &str) -> Result<Vec<AlignmentRecord>, Report> {
    Ok(
      read_many_fasta_str(fasta, &Alphabet::default())?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect(),
    )
  }

  fn node_name(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> String {
    names[&key].clone().expect("named node")
  }

  fn reconstruct_named_sparse(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    recon: &mut SparseReconstruction,
    impute: bool,
  ) -> Result<BTreeMap<String, Seq>, Report> {
    let mut out = BTreeMap::new();
    let SparseReconstruction {
      partition,
      node_states,
      edges,
      ..
    } = recon;
    let mut rng = rand::thread_rng();
    ancestral_reconstruction(graph, |node| {
      let seq = partition.reconstruct_node_sequence(
        node_states,
        &edges.forward,
        node,
        true,
        impute,
        SampleMode::Argmax,
        &mut rng,
      )?;
      out.insert(names[&node.key].clone().expect("named node"), seq);
      Some(())
    })?;
    Ok(out)
  }

  fn reconstruct_named_dense(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    recon: &mut DenseReconstruction,
    impute: bool,
  ) -> Result<BTreeMap<String, Seq>, Report> {
    let mut out = BTreeMap::new();
    let DenseReconstruction {
      partition, node_states, ..
    } = recon;
    let mut rng = rand::thread_rng();
    ancestral_reconstruction(graph, |node| {
      let seq = partition.reconstruct_node_sequence(node_states, node, true, impute, SampleMode::Argmax, &mut rng)?;
      out.insert(names[&node.key].clone().expect("named node"), seq);
      Some(())
    })?;
    Ok(out)
  }

  fn reconstruct_sparse(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    impute: bool,
  ) -> Result<BTreeMap<String, String>, Report> {
    let fitch = create_fitch_partition(
      graph,
      0,
      Alphabet::default(),
      &node_seq_inputs(graph, names, aln.to_vec()),
    )?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let recon = SparseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states);
    let (mut recon, _) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    Ok(to_strings(reconstruct_named_sparse(graph, names, &mut recon, impute)?))
  }

  fn reconstruct_dense(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    impute: bool,
  ) -> Result<BTreeMap<String, String>, Report> {
    let partition = PartitionMarginalDense::new(0, Alphabet::default(), get_common_length(aln)?);
    let node_states = partition.attach_sequences(graph, &node_seq_inputs(graph, names, aln.to_vec()))?;
    let recon = DenseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states);
    let (mut recon, _) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    Ok(to_strings(reconstruct_named_dense(graph, names, &mut recon, impute)?))
  }

  fn to_strings(seqs: BTreeMap<String, Seq>) -> BTreeMap<String, String> {
    seqs.into_iter().map(|(name, seq)| (name, seq.to_string())).collect()
  }
}
