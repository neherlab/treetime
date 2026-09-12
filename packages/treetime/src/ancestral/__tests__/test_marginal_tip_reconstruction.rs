#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::{
    ancestral_reconstruction_marginal, initialize_marginal, marginal_update, profile_branch_lengths,
  };
  use crate::ancestral::sample::SampleMode;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::traits::{PartitionBranchOps, PartitionMarginalOps};
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_primitives::Seq;

  /// C3: a tip that is Fitch-equal to its parent must keep its own observed nucleotide.
  ///
  /// `A=ACGT`, `B=GCGT` on `(A,B)`. The root MAP at position 0 is `G`, so the old parent-chained
  /// reconstruction emitted `A=GCGT` even though the reported mutation `G1A` was correct. The tip
  /// must echo the observed `ACGT`, and both backends must agree.
  #[test]
  fn test_marginal_tip_c3_preserves_observed_leaf() -> Result<(), Report> {
    let aln = parse_aln(indoc! {r#"
      >A
      ACGT
      >B
      GCGT
    "#})?;
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("(A:0.4,B:0.1)root:0.0;")?;
    let graph: Graph = graph;

    let sparse = reconstruct_sparse(&graph, &branch_lengths, &names, &aln, false)?;
    let dense = reconstruct_dense(&graph, &branch_lengths, &names, &aln, false)?;

    // Oracle: the observed input alignment. `A` retains its observed `ACGT`.
    assert_eq!("ACGT", sparse["A"]);
    assert_eq!("GCGT", sparse["B"]);
    assert_eq!(sparse, dense);
    Ok(())
  }

  /// F7 invariant: on every edge, applying the reported substitutions to the parent's reconstructed
  /// sequence reproduces the child's reconstructed sequence. The old tip corruption broke this on
  /// tip edges (dropping real parent->tip substitutions).
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
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("((A:0.1,B:0.1)AB:0.1,C:0.1)root:0.0;")?;
    let graph: Graph = graph;

    let alphabet = Alphabet::default();
    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln, &names)?;
    let mut partitions = [fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?];
    marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions)?;

    let seqs = reconstruct_named(&graph, &names, &mut partitions, false)?;

    let partition = &partitions[0];
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let parent = node_name(&names, edge.source());
      let child = node_name(&names, edge.target());
      let mut expected = seqs[&parent].clone();
      for sub in partition.edge_subs(&graph, edge.key())? {
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

  /// F8/C1: a tip `N` and a tip `R` echo the observed input without imputation and resolve to the
  /// parent-informed most likely state with imputation. Both backends must agree.
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
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("(A:0.1,(B:0.1,C:0.1)BC:0.1)root:0.0;")?;
    let graph: Graph = graph;

    // Without imputation the tip echoes its observed input (`N` and `R` preserved).
    let sparse_plain = reconstruct_sparse(&graph, &branch_lengths, &names, &aln, false)?;
    let dense_plain = reconstruct_dense(&graph, &branch_lengths, &names, &aln, false)?;
    assert_eq!("ANRT", sparse_plain["A"]);
    assert_eq!(sparse_plain, dense_plain);

    // With imputation, `N`->`C` and `R`->`G` (position 1 is `C` in every other tip; position 2 is
    // `G`, and `R={A,G}` resolves to `G`). Gaps would stay gaps; there are none here.
    let sparse_imputed = reconstruct_sparse(&graph, &branch_lengths, &names, &aln, true)?;
    let dense_imputed = reconstruct_dense(&graph, &branch_lengths, &names, &aln, true)?;
    assert_eq!("ACGT", sparse_imputed["A"]);
    assert_eq!(sparse_imputed, dense_imputed);
    Ok(())
  }

  /// Golden master vs TreeTime v0. With `--reconstruct-tip-states`, v0 imputes tip `A` of the
  /// alignment below to `ACGT` (captured from
  /// `./dev/docker/python treetime ancestral --reconstruct-tip-states`, which writes
  /// `>A\nACGT` to `ancestral_sequences.fasta`). v0 rejects two-taxon trees, so this uses three.
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
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("(A:0.1,(B:0.1,C:0.1)BC:0.1)root:0.0;")?;
    let graph: Graph = graph;

    let sparse = reconstruct_sparse(&graph, &branch_lengths, &names, &aln, true)?;
    let dense = reconstruct_dense(&graph, &branch_lengths, &names, &aln, true)?;

    // Oracle: TreeTime v0 Python reference output.
    assert_eq!("ACGT", sparse["A"]);
    assert_eq!("ACGT", dense["A"]);
    Ok(())
  }

  fn parse_aln(fasta: &str) -> Result<Vec<FastaRecord>, Report> {
    read_many_fasta_str(fasta, &Alphabet::default())
  }

  fn node_name(names: &BTreeMap<GraphNodeKey, Option<String>>, key: GraphNodeKey) -> String {
    names[&key].clone().expect("named node")
  }

  fn reconstruct_named<P>(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    partitions: &mut [P],
    impute: bool,
  ) -> Result<BTreeMap<String, Seq>, Report>
  where
    P: PartitionMarginalOps + crate::partition::traits::HasLogLh,
  {
    let mut out = BTreeMap::new();
    ancestral_reconstruction_marginal(
      graph,
      true,
      impute,
      partitions,
      SampleMode::Argmax,
      &mut rand::thread_rng(),
      |key, seq| {
        out.insert(names[&key].clone().expect("named node"), seq.clone());
        Ok(())
      },
    )?;
    Ok(out)
  }

  fn reconstruct_sparse(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[FastaRecord],
    impute: bool,
  ) -> Result<BTreeMap<String, String>, Report> {
    let fitch = create_fitch_partition(graph, 0, Alphabet::default(), aln, names)?;
    let mut partitions = [fitch.into_marginal_sparse(jc69(JC69Params::default())?, graph)?];
    marginal_update(graph, &profile_branch_lengths(branch_lengths), &mut partitions)?;
    Ok(to_strings(reconstruct_named(graph, names, &mut partitions, impute)?))
  }

  fn reconstruct_dense(
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[FastaRecord],
    impute: bool,
  ) -> Result<BTreeMap<String, String>, Report> {
    let partition = PartitionMarginalDense::new(
      0,
      jc69(JC69Params::default())?,
      Alphabet::default(),
      get_common_length(aln)?,
    );
    let mut partitions = [partition];
    initialize_marginal(
      graph,
      &profile_branch_lengths(branch_lengths),
      &mut partitions,
      aln,
      names,
    )?;
    Ok(to_strings(reconstruct_named(graph, names, &mut partitions, impute)?))
  }

  fn to_strings(seqs: BTreeMap<String, Seq>) -> BTreeMap<String, String> {
    seqs.into_iter().map(|(name, seq)| (name, seq.to_string())).collect()
  }
}
