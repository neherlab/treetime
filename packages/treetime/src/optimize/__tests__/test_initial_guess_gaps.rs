#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::{initialize_marginal, marginal_update, profile_branch_lengths};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::run_loop::optimize_partition_view;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use crate::partition::traits::PartitionBranchOps;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use treetime_graph::node::GraphNodeKey;

  use pretty_assertions::assert_eq;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::{NwkParse, nwk_read_str};

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGTACGTACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#},
      &alphabet,
    )
  }

  /// All four sequences share gaps at positions 4-7. Remaining 12 positions
  /// are identical to `gap_free_alignment`.
  fn gappy_alignment_shared() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGT----ACGTACGT
        >B
        ACGT----ACGTACGA
        >C
        ACGT----ACGTACGG
        >D
        ACGT----ACGTACGC
      "#},
      &alphabet,
    )
  }

  /// B alone has gaps at positions 4-7. Other sequences have ACGT there.
  fn gappy_alignment_one_leaf() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        ACGTACGTACGTACGT
        >B
        ACGT----ACGTACGA
        >C
        ACGTACGTACGTACGG
        >D
        ACGTACGTACGTACGC
      "#},
      &alphabet,
    )
  }

  fn setup_sparse(
    graph: &GraphAncestral,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[FastaRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<Vec<PartitionMarginalSparse>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let fitch = create_fitch_partition(graph, 0, alphabet, aln, names)?;
    let mut partitions = vec![fitch.into_marginal_sparse(jc69(JC69Params::default())?, graph)?];
    marginal_update(graph, &profile_branch_lengths(branch_lengths), &mut partitions)?.value();

    Ok(partitions)
  }

  fn setup_dense(
    graph: &GraphAncestral,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[FastaRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<Vec<PartitionMarginalDense>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let mut partitions = vec![PartitionMarginalDense::new(
      0,
      jc69(JC69Params::default())?,
      alphabet,
      get_common_length(aln)?,
    )];

    initialize_marginal(
      graph,
      &profile_branch_lengths(branch_lengths),
      &mut partitions,
      aln,
      names,
    )?
    .value();

    Ok(partitions)
  }

  fn get_branch_lengths(graph: &GraphAncestral, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|edge| branch_lengths[&edge.read_arc().key()].unwrap_or(0.0))
      .collect()
  }

  #[test]
  fn test_sparse_effective_length_no_gaps() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let graph: GraphAncestral = graph;
    let partitions = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = &partitions[0];
      let effective = p.edge_effective_length(&graph, edge_key)?;
      assert_eq!(16, effective);
    }

    Ok(())
  }

  #[test]
  fn test_dense_effective_length_no_gaps() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let graph: GraphAncestral = graph;
    let partitions = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = &partitions[0];
      let effective = p.edge_effective_length(&graph, edge_key)?;
      assert_eq!(16, effective);
    }

    Ok(())
  }

  #[test]
  fn test_sparse_effective_length_shared_gaps() -> Result<(), Report> {
    let aln = gappy_alignment_shared()?;
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let graph: GraphAncestral = graph;
    let partitions = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = &partitions[0];
      let effective = p.edge_effective_length(&graph, edge_key)?;
      // All nodes share gaps at positions 4-7, so effective = 16 - 4 = 12
      assert_eq!(12, effective);
    }

    Ok(())
  }

  #[test]
  fn test_dense_effective_length_shared_gaps() -> Result<(), Report> {
    let aln = gappy_alignment_shared()?;
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let graph: GraphAncestral = graph;
    let partitions = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = &partitions[0];
      let effective = p.edge_effective_length(&graph, edge_key)?;
      // All nodes share gaps at positions 4-7, so effective = 16 - 4 = 12
      assert_eq!(12, effective);
    }

    Ok(())
  }

  #[test]
  fn test_sparse_effective_length_one_leaf_gapped() -> Result<(), Report> {
    let aln = gappy_alignment_one_leaf()?;
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let graph: GraphAncestral = graph;
    let partitions = setup_sparse(&graph, &names, &aln, &branch_lengths)?;

    let mut found_reduced = false;
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let p = &partitions[0];
      let effective = p.edge_effective_length(&graph, edge_key)?;
      // At least one edge (B→AB) should have reduced effective length
      if effective < 16 {
        found_reduced = true;
      }
    }
    assert!(found_reduced);

    Ok(())
  }

  #[test]
  fn test_dense_edge_subs_excludes_gap_positions() -> Result<(), Report> {
    let aln = gappy_alignment_one_leaf()?;
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let graph: GraphAncestral = graph;
    let partitions = setup_dense(&graph, &names, &aln, &branch_lengths)?;

    let p = &partitions[0];
    for edge_ref in graph.get_edges() {
      let edge_key = edge_ref.read_arc().key();
      let subs = p.edge_subs(&graph, edge_key)?;
      // No substitution should involve a gap position (4-7)
      for sub in &subs {
        assert!(sub.pos() < 4 || sub.pos() >= 8);
      }
    }

    Ok(())
  }

  /// With shared gaps at positions 4-7 (which have identical nucleotides in
  /// the gap-free version), the initial guess should produce the same
  /// substitution rate per informative site. Branch lengths should be
  /// proportionally adjusted: subs/12 for gappy vs subs/16 for gap-free.
  #[test]
  fn test_initial_guess_sparse_gap_adjusted_rate() -> Result<(), Report> {
    let aln_clean = gap_free_alignment()?;
    let NwkParse {
      graph: graph_clean,
      names: graph_clean_names,
      branch_lengths: mut branch_lengths_clean,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let partitions_clean = setup_sparse(&graph_clean, &graph_clean_names, &aln_clean, &branch_lengths_clean)?;
    initial_guess_mixed(
      &graph_clean,
      &optimize_partition_view(&[], &partitions_clean),
      true,
      false,
      &mut branch_lengths_clean,
    )?;
    let bl_clean = get_branch_lengths(&graph_clean, &branch_lengths_clean);

    let aln_gappy = gappy_alignment_shared()?;
    let NwkParse {
      graph: graph_gappy,
      names: graph_gappy_names,
      branch_lengths: mut branch_lengths_gappy,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let partitions_gappy = setup_sparse(&graph_gappy, &graph_gappy_names, &aln_gappy, &branch_lengths_gappy)?;
    initial_guess_mixed(
      &graph_gappy,
      &optimize_partition_view(&[], &partitions_gappy),
      true,
      false,
      &mut branch_lengths_gappy,
    )?;
    let bl_gappy = get_branch_lengths(&graph_gappy, &branch_lengths_gappy);

    // With 4 shared gap positions out of 16, the effective length is 12.
    // Substitutions at non-gap positions are the same, so the per-site rate
    // is higher by factor 16/12 = 4/3.
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
