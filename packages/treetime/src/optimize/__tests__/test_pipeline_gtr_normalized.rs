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
  use crate::test_utils::{TestDataset, leaf_seq_inputs, sequence_source};
  use eyre::Report;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_optimize_pipeline_infer_gtr_mu_normalized() -> Result<(), Report> {
    let alphabet = Alphabet::default();
    let nwk_parsed = nwk_read(TestDataset::Flu20.tree().as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let sequences = leaf_seq_inputs(&graph, &names, TestDataset::Flu20.alignment()?);

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
      names,
      sources: vec![sequence_source(alphabet, sequences)?],
      branch_lengths,
    };

    let output = run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let [result] = output.results.as_slice() else {
      panic!("one alignment gives one result");
    };
    let gtr = result.gtr.as_ref().expect("optimize fits a model");
    pretty_assert_ulps_eq!(gtr.mu(), 1.0, max_ulps = 4);

    Ok(())
  }
}
