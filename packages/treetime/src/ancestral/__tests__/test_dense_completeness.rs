#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::partition::traits::PartitionBranchOps;
  use crate::partition::traits::PartitionOptimizeOps;
  use crate::partition::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  fn setup_dense_with_unknowns() -> Result<(GraphAncestral, Arc<RwLock<PartitionMarginalDense>>), Report> {
    let newick = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let fasta = "
>A
ACGTNNACGT
>B
ACGTACNNGT
>C
ACGTACGTNN
>D
NNGTACGTAC
";
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let length = get_common_length(&aln)?;

    let partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    initialize_marginal(&graph, std::slice::from_ref(&partition), &aln)?;
    Ok((graph, partition))
  }

  fn setup_sparse_with_unknowns() -> Result<(GraphAncestral, Arc<RwLock<PartitionMarginalSparse>>), Report> {
    let newick = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let fasta = "
>A
ACGTNNACGT
>B
ACGTACNNGT
>C
ACGTACGTNN
>D
NNGTACGTAC
";
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let length = get_common_length(&aln)?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
    let partition = Arc::new(RwLock::new(
      fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?,
    ));
    update_marginal(&graph, std::slice::from_ref(&partition))?;
    Ok((graph, partition))
  }

  #[test]
  fn test_dense_completeness_non_char_tracked_on_leaves() -> Result<(), Report> {
    let (_, partition) = setup_dense_with_unknowns()?;
    let p = partition.read_arc();

    // Find leaf A's node (has "NN" at positions 4-5)
    let leaf_a = p.nodes.values().find(|n| n.seq.unknown.contains(&(4, 6)));
    assert!(leaf_a.is_some(), "Leaf A should have unknown range (4,6)");

    let leaf_a = leaf_a.unwrap();
    assert!(!leaf_a.seq.non_char.is_empty(), "Leaf A should have non_char ranges");
    assert!(
      leaf_a.seq.non_char.contains(&(4, 6)),
      "Leaf A non_char should contain (4,6)"
    );
    Ok(())
  }

  #[test]
  fn test_dense_completeness_effective_length_subtracts_unknowns() -> Result<(), Report> {
    let (graph, partition) = setup_dense_with_unknowns()?;
    let p = partition.read_arc();

    // Every leaf has 2 N positions. Leaf edges have the leaf's unknowns in
    // their non_char, reducing effective length below 10. Internal edges may
    // still have effective length 10 if children's unknown positions don't overlap.
    let mut any_reduced = false;
    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      let eff_len = p.edge_effective_length(&graph, edge_key)?;
      assert!(
        eff_len <= 10,
        "Effective length {eff_len} should not exceed total length 10"
      );
      if eff_len < 10 {
        any_reduced = true;
      }
    }
    assert!(
      any_reduced,
      "At least one edge should have reduced effective length due to unknown positions"
    );
    Ok(())
  }

  #[test]
  fn test_dense_sparse_effective_length_agreement() -> Result<(), Report> {
    let (graph_d, partition_d) = setup_dense_with_unknowns()?;
    let (graph_s, partition_s) = setup_sparse_with_unknowns()?;

    let pd = partition_d.read_arc();
    let ps = partition_s.read_arc();

    // Both graphs have same topology, edges in same order
    let dense_edges = graph_d.get_edges();
    let sparse_edges = graph_s.get_edges();
    assert_eq!(dense_edges.len(), sparse_edges.len());

    for (de, se) in dense_edges.iter().zip(sparse_edges.iter()) {
      let dk = de.read_arc().key();
      let sk = se.read_arc().key();

      let dense_eff = pd.edge_effective_length(&graph_d, dk)?;
      let sparse_eff = ps.edge_effective_length(&graph_s, sk)?;

      assert_eq!(
        dense_eff, sparse_eff,
        "Edge effective length mismatch: dense={dense_eff}, sparse={sparse_eff} for edge {dk:?}/{sk:?}"
      );
    }

    Ok(())
  }

  fn setup_dense_with_gaps() -> Result<(GraphAncestral, Arc<RwLock<PartitionMarginalDense>>), Report> {
    let newick = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let fasta = "
>A
ACGT--ACGT
>B
ACGTACACGT
>C
AC--ACGTAC
>D
ACGTACGTAC
";
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let length = get_common_length(&aln)?;

    let partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length,
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    initialize_marginal(&graph, std::slice::from_ref(&partition), &aln)?;
    Ok((graph, partition))
  }

  #[test]
  fn test_dense_completeness_indels_populated() -> Result<(), Report> {
    let (graph, partition) = setup_dense_with_gaps()?;
    let p = partition.read_arc();

    // Alignment: A=ACGT--ACGT, B=ACGTACACGT, C=AC--ACGTAC, D=ACGTACGTAC
    // A has gap at (4,6), C has gap at (2,4). B and D have no gaps.
    // Indels should appear on edges connecting to A and C.
    let total_indels: usize = p.edges.values().map(|e| e.indels.len()).sum();
    assert!(
      total_indels >= 2,
      "Expected at least 2 indels (one for A's gap, one for C's gap), got {total_indels}"
    );

    // Verify indel directions exist (at least one deletion)
    let has_deletion = p.edges.values().any(|e| e.indels.iter().any(|i| i.deletion));
    assert!(
      has_deletion,
      "At least one deletion should be detected from gap-bearing leaves"
    );
    Ok(())
  }

  #[test]
  fn test_dense_completeness_edge_indel_count_nonzero() -> Result<(), Report> {
    let (graph, partition) = setup_dense_with_gaps()?;
    let p = partition.read_arc();

    let total: usize = graph
      .get_edges()
      .iter()
      .map(|e| p.edge_indel_count(e.read_arc().key()))
      .sum();

    assert!(
      total >= 2,
      "edge_indel_count total should be at least 2 for gaps at (4,6) and (2,4), got {total}"
    );
    Ok(())
  }

  fn setup_sparse_with_gaps() -> Result<(GraphAncestral, Arc<RwLock<PartitionMarginalSparse>>), Report> {
    let newick = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let fasta = "
>A
ACGT--ACGT
>B
ACGTACACGT
>C
AC--ACGTAC
>D
ACGTACGTAC
";
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let length = get_common_length(&aln)?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln)?;
    let partition = Arc::new(RwLock::new(
      fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?,
    ));
    update_marginal(&graph, std::slice::from_ref(&partition))?;
    Ok((graph, partition))
  }

  #[test]
  fn test_dense_sparse_indel_count_agreement() -> Result<(), Report> {
    let (graph_d, partition_d) = setup_dense_with_gaps()?;
    let (graph_s, partition_s) = setup_sparse_with_gaps()?;

    let pd = partition_d.read_arc();
    let ps = partition_s.read_arc();

    let dense_total: usize = graph_d
      .get_edges()
      .iter()
      .map(|e| pd.edge_indel_count(e.read_arc().key()))
      .sum();

    let sparse_total: usize = graph_s
      .get_edges()
      .iter()
      .map(|e| ps.edge_indel_count(e.read_arc().key()))
      .sum();

    assert_eq!(
      dense_total, sparse_total,
      "Total indel count should agree between dense ({dense_total}) and sparse ({sparse_total})"
    );
    Ok(())
  }
}
