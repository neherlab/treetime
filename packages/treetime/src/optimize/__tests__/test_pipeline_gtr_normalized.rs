#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
  use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
  use crate::partition::traits::HasGtr;
  use crate::payload::ancestral::GraphAncestral;
  use crate::progress::NoopProgress;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use std::path::Path;
  use treetime_io::fasta::read_many_fasta;
  use treetime_io::nwk::nwk_read_file;

  // `optimize --gtr=infer` must serialize the GTR after rate normalization.
  // `normalize_partition_rates` rescales `mu` so the average substitution rate
  // is 1 and absorbs the rate into branch lengths; for a single partition this
  // drives `mu` to exactly 1.0. Reading a creation-time snapshot instead would
  // report the raw inferred `mu` (essentially never 1.0) while shipping
  // rate-scaled branch lengths, an internally inconsistent model. This test
  // fails on that stale snapshot.
  #[test]
  fn test_optimize_pipeline_infer_gtr_mu_normalized() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(Path::parent)
      .expect("workspace root");

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = workspace_root.join("data/flu/h3n2/20/aln.fasta.xz");

    let graph: GraphAncestral = nwk_read_file(&tree_path)?;
    let sequences = read_many_fasta(&[aln_path.to_str().expect("utf-8 path")], &alphabet)?;

    let params = OptimizeParams {
      model: GtrModelName::Infer,
      dense: Some(false),
      max_iter: 2,
      dp: 0.1,
      damping: 0.0,
      opt_method: BranchOptMethod::default(),
      initial_guess: InitialGuessMode::default(),
      no_indels: false,
      reroot_spec: None,
    };

    let input = OptimizeInput {
      graph,
      alphabet,
      sequences,
    };

    let output = run(&params, input, &NoopProgress)?;

    // Normalization contract: average rate 1 => mu == 1.0 for a single partition.
    assert_ulps_eq!(output.gtr.mu, 1.0, max_ulps = 4);

    // Single source of truth: the serialized GTR is the partition's live model,
    // not an independent snapshot.
    let partition_mu = output.sparse_partitions[0].read_arc().gtr().mu;
    assert_ulps_eq!(output.gtr.mu, partition_mu, max_ulps = 4);

    Ok(())
  }
}
