#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use crate::seq::indel::InDel;
  use eyre::Report;
  use treetime_graph::graph::Graph;

  use pretty_assertions::assert_eq;

  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::{NwkParse, nwk_read_str};

  fn setup_dense_with_unknowns() -> Result<(Graph, DenseReconstruction), Report> {
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
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(newick)?;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let length = get_common_length(&aln)?;

    let partition = PartitionMarginalDense::new(0, jc69(JC69Params::default())?, alphabet, length);
    let node_states = partition.attach_sequences(&graph, &aln, &names)?;
    let recon = DenseReconstruction::seeded(partition, node_states);
    let (recon, _) = recon.marginal_update(&graph, &profile_branch_lengths(&branch_lengths))?;
    Ok((graph, recon))
  }

  fn setup_sparse_with_unknowns() -> Result<(Graph, SparseReconstruction), Report> {
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
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(newick)?;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let length = get_common_length(&aln)?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln, &names)?;
    let (partition, node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let recon = SparseReconstruction::seeded(partition, node_states);
    let (recon, _) = recon.marginal_update(&graph, &profile_branch_lengths(&branch_lengths))?;
    Ok((graph, recon))
  }

  #[test]
  fn test_dense_completeness_non_char_tracked_on_leaves() -> Result<(), Report> {
    let (_, recon) = setup_dense_with_unknowns()?;

    // Find leaf A's node (has "NN" at positions 4-5)
    let leaf_a = recon.node_states.values().find(|n| n.seq.unknown.contains(&(4, 6)));
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
    let (graph, recon) = setup_dense_with_unknowns()?;

    // Every leaf has 2 N positions. Leaf edges have the leaf's unknowns in
    // their non_char, reducing effective length below 10. Internal edges may
    // still have effective length 10 if children's unknown positions don't overlap.
    let mut any_reduced = false;
    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      let eff_len = recon
        .partition
        .edge_effective_length(&recon.node_states, &graph, edge_key)?;
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
    let (graph_d, recon_d) = setup_dense_with_unknowns()?;
    let (graph_s, recon_s) = setup_sparse_with_unknowns()?;

    // Both graphs have same topology, edges in same order
    let dense_edges = graph_d.get_edges();
    let sparse_edges = graph_s.get_edges();
    assert_eq!(dense_edges.len(), sparse_edges.len());

    for (de, se) in dense_edges.iter().zip(sparse_edges.iter()) {
      let dk = de.read_arc().key();
      let sk = se.read_arc().key();

      let dense_eff = recon_d
        .partition
        .edge_effective_length(&recon_d.node_states, &graph_d, dk)?;
      let sparse_eff = recon_s.partition.edge_effective_length(&graph_s, sk)?;

      assert_eq!(
        dense_eff, sparse_eff,
        "Edge effective length mismatch: dense={dense_eff}, sparse={sparse_eff} for edge {dk:?}/{sk:?}"
      );
    }

    Ok(())
  }

  fn setup_dense_with_gaps() -> Result<(Graph, DenseReconstruction), Report> {
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
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(newick)?;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let length = get_common_length(&aln)?;

    let partition = PartitionMarginalDense::new(0, jc69(JC69Params::default())?, alphabet, length);
    let node_states = partition.attach_sequences(&graph, &aln, &names)?;
    let recon = DenseReconstruction::seeded(partition, node_states);
    let (recon, _) = recon.marginal_update(&graph, &profile_branch_lengths(&branch_lengths))?;
    Ok((graph, recon))
  }

  #[test]
  fn test_dense_completeness_indels_populated() -> Result<(), Report> {
    let (_graph, recon) = setup_dense_with_gaps()?;

    // Alignment: A=ACGT--ACGT, B=ACGTACACGT, C=AC--ACGTAC, D=ACGTACGTAC
    // A has gap at (4,6), C has gap at (2,4). B and D have no gaps.
    // Indels should appear on edges connecting to A and C.
    let total_indels: usize = recon.edges.estimates.values().map(|e| e.indels.len()).sum();
    assert!(
      total_indels >= 2,
      "Expected at least 2 indels (one for A's gap, one for C's gap), got {total_indels}"
    );

    // Verify indel directions exist (at least one deletion)
    let has_deletion = recon
      .edges
      .estimates
      .values()
      .any(|e| e.indels.iter().any(InDel::is_deletion));
    assert!(
      has_deletion,
      "At least one deletion should be detected from gap-bearing leaves"
    );
    Ok(())
  }

  #[test]
  fn test_dense_completeness_edge_indel_count_nonzero() -> Result<(), Report> {
    let (graph, recon) = setup_dense_with_gaps()?;

    let total: usize = graph
      .get_edges()
      .iter()
      .map(|e| {
        recon
          .partition
          .edge_indel_count(&recon.edges.estimates, e.read_arc().key())
      })
      .sum();

    assert!(
      total >= 2,
      "edge_indel_count total should be at least 2 for gaps at (4,6) and (2,4), got {total}"
    );
    Ok(())
  }

  fn setup_sparse_with_gaps() -> Result<(Graph, SparseReconstruction), Report> {
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
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(newick)?;
    let graph: Graph = graph;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let length = get_common_length(&aln)?;

    let fitch = create_fitch_partition(&graph, 0, alphabet, &aln, &names)?;
    let (partition, node_states) = fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?;
    let recon = SparseReconstruction::seeded(partition, node_states);
    let (recon, _) = recon.marginal_update(&graph, &profile_branch_lengths(&branch_lengths))?;
    Ok((graph, recon))
  }

  #[test]
  fn test_dense_sparse_indel_count_agreement() -> Result<(), Report> {
    let (graph_d, recon_d) = setup_dense_with_gaps()?;
    let (graph_s, recon_s) = setup_sparse_with_gaps()?;

    let dense_total: usize = graph_d
      .get_edges()
      .iter()
      .map(|e| {
        recon_d
          .partition
          .edge_indel_count(&recon_d.edges.estimates, e.read_arc().key())
      })
      .sum();

    let sparse_total: usize = graph_s
      .get_edges()
      .iter()
      .map(|e| recon_s.partition.edge_indel_count(e.read_arc().key()))
      .sum();

    assert_eq!(
      dense_total, sparse_total,
      "Total indel count should agree between dense ({dense_total}) and sparse ({sparse_total})"
    );
    Ok(())
  }
}
