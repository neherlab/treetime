#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::commands::optimize::args::BranchOptMethod;
  use crate::commands::optimize::optimize_unified::initial_guess_mixed;
  use crate::commands::optimize::run::{collect_optimize_partitions, run_optimize_loop};
  use crate::gtr::get_gtr::{JC69Params, get_gtr_sparse, jc69};
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::path::Path;
  use std::sync::Arc;
  use treetime_io::fasta::read_many_fasta;
  use treetime_io::nwk::nwk_read_file;
  use treetime_primitives::seq;

  /// Regression test: sparse optimize loop converges on sc2/2844 (dataset with indels).
  ///
  /// This dataset exhibits a persistent likelihood 2-cycle without the three-condition convergence check.
  /// The loop alternated between ~-143156 and ~-143157 from iteration ~15 onward,
  /// exhausting all max_iter iterations. With the three-condition convergence check,
  /// damping floor, restored indel rate estimation, and v0-aligned defaults, the loop
  /// must stop via one of the three conditions before exhausting max_iter=50.
  #[test]
  fn test_convergence_sc2_sparse_converges_on_sc2_2844() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/sc2/2844/tree.nwk");
    let aln_path = workspace_root.join("data/sc2/2844/aln.fasta.xz");
    let aln = read_many_fasta(&[aln_path.to_str().unwrap()], &alphabet)?;
    let mut graph = nwk_read_file(&tree_path)?;

    let sparse_partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }))];

    compress_sequences(&graph, &sparse_partitions, &aln)?;
    for p in &sparse_partitions {
      p.write_arc().extract_root_sequence(&graph)?;
    }
    update_marginal(&graph, &sparse_partitions)?;

    // Infer GTR from data (matches production flow)
    for partition in &sparse_partitions {
      let gtr = get_gtr_sparse(&crate::gtr::get_gtr::GtrModelName::JC69, partition, &graph)?;
      partition.write_arc().gtr = gtr;
    }

    let dense_partitions = vec![];
    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);
    initial_guess_mixed(&graph, &mixed_partitions, true)?;

    let max_iter = 50;
    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      max_iter,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
    )?;

    assert!(
      result.stopped_at.is_some(),
      "Sparse optimize loop on sc2/2844 did not converge within {max_iter} iterations. LH history (last 5): {:?}",
      result.lh_history.iter().rev().take(5).collect::<Vec<_>>()
    );

    Ok(())
  }

  /// Regression test: flu/h3n2/20 converges normally (indel-free dataset).
  ///
  /// Verifies the fix does not regress convergence on a dataset where the
  /// variable position set is stable (no oscillation).
  #[test]
  fn test_convergence_sc2_flu_h3n2_20_converges() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .unwrap();

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = workspace_root.join("data/flu/h3n2/20/aln.fasta.xz");
    let aln = read_many_fasta(&[aln_path.to_str().unwrap()], &alphabet)?;
    let mut graph = nwk_read_file(&tree_path)?;

    let sparse_partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }))];

    compress_sequences(&graph, &sparse_partitions, &aln)?;
    for p in &sparse_partitions {
      p.write_arc().extract_root_sequence(&graph)?;
    }
    update_marginal(&graph, &sparse_partitions)?;

    for partition in &sparse_partitions {
      let gtr = get_gtr_sparse(&crate::gtr::get_gtr::GtrModelName::JC69, partition, &graph)?;
      partition.write_arc().gtr = gtr;
    }

    let dense_partitions = vec![];
    let mixed_partitions = collect_optimize_partitions(&dense_partitions, &sparse_partitions);
    initial_guess_mixed(&graph, &mixed_partitions, true)?;

    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      10,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
    )?;

    assert!(
      result.stopped_at.is_some(),
      "flu/h3n2/20 did not converge within 10 iterations"
    );
    Ok(())
  }
}
