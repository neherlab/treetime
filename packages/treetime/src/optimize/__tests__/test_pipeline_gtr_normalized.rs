#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::cancel::NoopCancel;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::OptimizeLoopSettings;
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
  use crate::optimize::partition_rates::PartitionRates;
  use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
  use crate::progress::NoopProgress;
  use crate::test_utils::{leaf_seq_inputs, sequence_source};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use std::path::Path;
  use treetime_io::fasta::fasta_read_file;
  use treetime_io::nwk::nwk_read_file;
  use treetime_primitives::AlignmentRecord;

  #[test]
  fn test_optimize_pipeline_infer_gtr_mu_normalized() -> Result<(), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(Path::parent)
      .expect("workspace root");

    let alphabet = Alphabet::default();
    let tree_path = workspace_root.join("data/flu/h3n2/20/tree.nwk");
    let aln_path = workspace_root.join("data/flu/h3n2/20/aln.fasta.xz");

    let nwk_parsed = nwk_read_file(&tree_path)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let sequences: Vec<AlignmentRecord> = fasta_read_file(&aln_path, &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    let sequences = leaf_seq_inputs(&graph, &names, sequences);

    let params = OptimizeParams {
      model: GtrModelName::Infer,
      dense: Some(false),
      initial_guess: InitialGuessMode::default(),
      reroot_spec: None,
      loop_settings: OptimizeLoopSettings {
        max_iter: 2,
        dp: 0.1,
        damping: 0.0,
        opt_method: BranchOptMethod::default(),
        indel_policy: IndelPolicy::Include,
        topology_ops: TopologyOps::default(),
        partition_rates: PartitionRates::Estimate,
      },
    };
    let input = OptimizeInput {
      graph,
      sources: vec![sequence_source(alphabet, sequences)?],
      branch_lengths,
    };

    let output = run(&params, input, &names, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let [result] = output.results.as_slice() else {
      panic!("one alignment gives one result");
    };
    let gtr = result.gtr.as_ref().expect("optimize fits a model");
    assert_ulps_eq!(gtr.mu(), 1.0, max_ulps = 4);

    Ok(())
  }
}
