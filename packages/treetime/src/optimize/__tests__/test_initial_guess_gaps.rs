#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::guess_initial_branch_lengths;
  use crate::optimize::gather::InitialGuessEvidence;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::ExistingBranchLengths;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, SparseReconstruction};
  use crate::partition::set::PartitionSet;
  use crate::test_utils::SingleReconstruction;
  use crate::test_utils::leaf_seq_inputs;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::{FastaRecord, fasta_read};
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    fasta_read(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGTACGTACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#}
      .as_bytes(),
      &alphabet,
    )
  }

  fn gappy_alignment_shared() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    fasta_read(
      indoc! {r#"
        >A
        ACGT----ACGTACGT
        >B
        ACGT----ACGTACGA
        >C
        ACGT----ACGTACGG
        >D
        ACGT----ACGTACGC
      "#}
      .as_bytes(),
      &alphabet,
    )
  }

  fn gappy_alignment_one_leaf() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    fasta_read(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGT----ACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#}
      .as_bytes(),
      &alphabet,
    )
  }

  fn setup_sparse(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<PartitionSet, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = create_fitch_partition(graph, alphabet, leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let (partition, node_states) = fitch.into_marginal_sparse(graph)?;
    let reconstruction = PartitionSet::Sparse(vec![SparseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
      node_states,
    )]);
    let (reconstruction, _) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;

    Ok(reconstruction)
  }

  fn setup_dense(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<PartitionSet, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let partition = PartitionMarginalDense::new(alphabet, graph, &leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let reconstruction = PartitionSet::Dense(vec![DenseReconstruction::seeded(
      partition,
      jc69(JC69Params::default())?,
    )]);

    let (reconstruction, _) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;

    Ok(reconstruction)
  }

  fn get_branch_lengths(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    graph
      .get_edges()
      .map(|edge| branch_lengths[&edge.key()].unwrap_or(0.0))
      .collect()
  }

  #[test]
  fn test_sparse_effective_length_no_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gap_free_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = reconstruction.single().edge_effective_length(&graph, edge_key)?;
      assert_eq!(16, effective);
    }

    Ok(())
  }

  #[test]
  fn test_dense_effective_length_no_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gap_free_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = reconstruction.single().edge_effective_length(&graph, edge_key)?;
      assert_eq!(16, effective);
    }

    Ok(())
  }

  #[test]
  fn test_sparse_effective_length_shared_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gappy_alignment_shared()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = reconstruction.single().edge_effective_length(&graph, edge_key)?;
      assert_eq!(12, effective);
    }

    Ok(())
  }

  #[test]
  fn test_dense_effective_length_shared_gaps() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gappy_alignment_shared()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = reconstruction.single().edge_effective_length(&graph, edge_key)?;
      assert_eq!(12, effective);
    }

    Ok(())
  }

  #[test]
  fn test_sparse_effective_length_one_leaf_gapped() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gappy_alignment_one_leaf()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    let mut found_reduced = false;
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let effective = reconstruction.single().edge_effective_length(&graph, edge_key)?;
      if effective < 16 {
        found_reduced = true;
      }
    }
    assert!(found_reduced);

    Ok(())
  }

  #[test]
  fn test_dense_edge_subs_excludes_gap_positions() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = gappy_alignment_one_leaf()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.key();
      let subs = reconstruction.single().edge_subs(&graph, edge_key)?;
      for sub in &subs {
        assert!(sub.pos() < 4 || sub.pos() >= 8);
      }
    }

    Ok(())
  }

  #[test]
  fn test_initial_guess_sparse_gap_adjusted_rate() -> Result<(), Report> {
    let aln_clean: Vec<AlignmentRecord> = gap_free_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph_clean_names = nwk_parsed.names();
    let graph_clean = nwk_parsed.graph;
    let mut branch_lengths_clean = nwk_parsed.branch_lengths;
    let reconstruction_clean = setup_sparse(&graph_clean, &graph_clean_names, &aln_clean, &branch_lengths_clean)?;
    {
      let total_length = reconstruction_clean.total_length();
      let indels = IndelPolicy::Include.edge_indels(&graph_clean, &reconstruction_clean, &branch_lengths_clean);
      let evidence = InitialGuessEvidence::gather(&graph_clean, &reconstruction_clean)?;
      guess_initial_branch_lengths(
        &graph_clean,
        &indels,
        &evidence,
        ExistingBranchLengths::Overwrite,
        &mut branch_lengths_clean,
      )?;
    }
    let bl_clean = get_branch_lengths(&graph_clean, &branch_lengths_clean);

    let aln_gappy: Vec<AlignmentRecord> = gappy_alignment_shared()?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph_gappy_names = nwk_parsed.names();
    let graph_gappy = nwk_parsed.graph;
    let mut branch_lengths_gappy = nwk_parsed.branch_lengths;
    let reconstruction_gappy = setup_sparse(&graph_gappy, &graph_gappy_names, &aln_gappy, &branch_lengths_gappy)?;
    {
      let indels = IndelPolicy::Include.edge_indels(&graph_gappy, &reconstruction_gappy, &branch_lengths_gappy);
      let evidence = InitialGuessEvidence::gather(&graph_gappy, &reconstruction_gappy)?;
      guess_initial_branch_lengths(
        &graph_gappy,
        &indels,
        &evidence,
        ExistingBranchLengths::Overwrite,
        &mut branch_lengths_gappy,
      )?;
    }
    let bl_gappy = get_branch_lengths(&graph_gappy, &branch_lengths_gappy);

    let ratio = 16.0 / 12.0;
    for (clean, gappy) in bl_clean.iter().zip(bl_gappy.iter()) {
      if *clean > 0.0 {
        assert_abs_diff_eq!(gappy / clean, ratio, epsilon = 1e-10);
      } else {
        assert_abs_diff_eq!(*gappy, 0.0, epsilon = 1e-10);
      }
    }

    Ok(())
  }
}
