#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction, SparseReconstruction};
  use crate::partition::marginal::sample::SampleMode;
  use crate::test_utils::leaf_seq_inputs;
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
  fn test_marginal_tip_c3_preserves_observed_leaf() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >A
      ACGT
      >B
      GCGT
    "#})?;
    let nwk_parsed = nwk_read(b"(A:0.4,B:0.1)root:0.0;".as_slice())?;
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
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.1)AB:0.1,C:0.1)root:0.0;".as_slice())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let alphabet = Alphabet::default();
    let fitch = create_fitch_partition(&graph, alphabet, leaf_seq_inputs(&graph, &names, aln))?;
    let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
    let recon = MarginalReconstruction::Sparse(SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    ));
    let (recon, _) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let seqs = reconstruct_named(&graph, &names, &recon, false)?;

    for edge in graph.get_edges() {
      let parent = node_name(&names, edge.source());
      let child = node_name(&names, edge.target());
      let mut expected = seqs[&parent].clone();
      for sub in recon.edge_subs(&graph, edge.key())? {
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
    let nwk_parsed = nwk_read(b"(A:0.1,(B:0.1,C:0.1)BC:0.1)root:0.0;".as_slice())?;
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
  #[cfg_attr(
    dylint_lib = "custom",
    expect(versioned_name, reason = "v0 names the reference implementation, not a revision")
  )]
  fn test_gm_marginal_tip_impute_matches_v0() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >A
      ANRT
      >B
      ACGT
      >C
      ACGT
    "#})?;
    let nwk_parsed = nwk_read(b"(A:0.1,(B:0.1,C:0.1)BC:0.1)root:0.0;".as_slice())?;
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
      fasta_read(fasta.as_bytes(), &Alphabet::default())?
        .into_iter()
        .map(AlignmentRecord::from)
        .collect(),
    )
  }

  fn node_name(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> String {
    names[&key].clone().expect("named node")
  }

  fn reconstruct_named(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    recon: &MarginalReconstruction,
    impute: bool,
  ) -> Result<BTreeMap<String, Seq>, Report> {
    let sampled = recon.sample_sequences(graph, SampleMode::Argmax, &mut get_random_number_generator(0))?;
    emitted_sequences_by_name(names, &node_keys(graph), &sampled, |key| {
      recon.node_sequence(graph, impute, key)
    })
  }

  fn reconstruct_sparse(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    impute: bool,
  ) -> Result<BTreeMap<String, String>, Report> {
    let fitch = create_fitch_partition(graph, Alphabet::default(), leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let recon = MarginalReconstruction::Sparse(SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    ));
    let (recon, _) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    Ok(to_strings(reconstruct_named(graph, names, &recon, impute)?))
  }

  fn reconstruct_dense(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    impute: bool,
  ) -> Result<BTreeMap<String, String>, Report> {
    let partition =
      PartitionMarginalDense::new(Alphabet::default(), graph, &leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let recon = MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, jc69(JC69Params::default())?));
    let (recon, _) = recon.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    Ok(to_strings(reconstruct_named(graph, names, &recon, impute)?))
  }

  fn to_strings(seqs: BTreeMap<String, Seq>) -> BTreeMap<String, String> {
    seqs.into_iter().map(|(name, seq)| (name, seq.to_string())).collect()
  }
}
