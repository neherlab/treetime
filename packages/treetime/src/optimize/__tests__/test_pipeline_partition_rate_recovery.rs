#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{FILE_COUNT, REPLICATE_COUNT, SIMULATED_RATES, estimated_rates};
  use itertools::{Itertools, izip};
  use ndarray::{Array1, Axis};
  use pretty_assertions::assert_eq;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::make_internal_report;

  #[test]
  fn test_optimize_estimated_partition_rates_recover_the_simulated_relative_rates() -> Result<(), Report> {
    let simulated = Array1::from_iter(SIMULATED_RATES);
    let targets = &simulated
      / simulated
        .mean()
        .ok_or_else(|| make_internal_report!("No simulated rates"))?;

    let estimates = estimated_rates()?;
    let means = estimates
      .mean_axis(Axis(0))
      .ok_or_else(|| make_internal_report!("No replicates"))?;
    let standard_errors = estimates.std_axis(Axis(0), 1.0) / REPLICATE_COUNT.to_f64_exact().sqrt();

    let within_three_standard_errors = izip!(&means, &standard_errors, &targets)
      .map(|(mean, standard_error, target)| (mean - target).abs() <= 3.0 * standard_error)
      .collect_vec();
    assert_eq!(
      vec![true; FILE_COUNT],
      within_three_standard_errors,
      "means {means}, standard errors {standard_errors}, targets {targets}"
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
    use crate::progress::NoopProgress;
    use crate::seq::alignment::SequenceSource;
    use crate::test_utils::{TestDataset, columns, fasta_xz_read, leaf_seq_inputs, sequence_source};
    use eyre::Report;
    use ndarray::Array2;
    use rayon::ThreadPoolBuilder;
    use rayon::iter::{IntoParallelIterator, ParallelIterator};
    use std::path::PathBuf;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::AlignmentRecord;
    use treetime_utils::make_internal_report;

    pub(super) const SIMULATED_RATES: [f64; 3] = [0.3, 1.0, 2.0];
    pub(super) const FILE_COUNT: usize = SIMULATED_RATES.len();
    pub(super) const REPLICATE_COUNT: usize = 20;
    const SITES_PER_FILE: usize = 3000;
    const THREAD_COUNT: usize = 8;
    const ALIGNMENTS_XZ: &[u8] = include_bytes!("__fixtures__/rate_recovery_alignments.fasta.xz");

    pub(super) fn estimated_rates() -> Result<Array2<f64>, Report> {
      let records = fasta_xz_read(ALIGNMENTS_XZ)?;
      let rates = ThreadPoolBuilder::new()
        .num_threads(THREAD_COUNT)
        .build()?
        .install(|| {
          (0..REPLICATE_COUNT)
            .into_par_iter()
            .map(|replicate| replicate_rates(&records, replicate))
            .collect::<Result<Vec<_>, Report>>()
        })?;
      Ok(Array2::from_shape_vec((REPLICATE_COUNT, FILE_COUNT), rates.concat())?)
    }

    fn replicate_rates(records: &[AlignmentRecord], replicate: usize) -> Result<Vec<f64>, Report> {
      let parse = nwk_read(TestDataset::Zika86.tree().as_bytes())?;
      let names = parse.names.clone();
      let sources = (0..FILE_COUNT)
        .map(|file| {
          let start = (replicate * FILE_COUNT + file) * SITES_PER_FILE;
          let source = sequence_source(
            Alphabet::default(),
            leaf_seq_inputs(&parse.graph, &names, columns(records, start..start + SITES_PER_FILE)),
          )?;
          Ok(SequenceSource {
            path: PathBuf::from(format!("{}.fasta", file + 1)),
            ..source
          })
        })
        .collect::<Result<Vec<_>, Report>>()?;
      let params = OptimizeParams {
        model: GtrModelName::JC69,
        dense: Some(false),
        initial_guess: InitialGuessMode::Auto,
        reroot_spec: None,
        loop_settings: OptimizeLoopSettings {
          max_iter: 10,
          dp: 0.1,
          damping: 0.0,
          opt_method: BranchOptMethod::Brent,
          indel_policy: IndelPolicy::Include,
          topology_ops: TopologyOps::default(),
          partition_rates: PartitionRates::Estimate,
        },
      };
      let input = OptimizeInput {
        graph: parse.graph,
        names,
        sources,
        branch_lengths: parse.branch_lengths,
      };
      let output = run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress)?;
      output
        .results
        .iter()
        .map(|result| {
          result
            .gtr
            .as_ref()
            .map(|gtr| gtr.mu())
            .ok_or_else(|| make_internal_report!("Partition {} has no model", result.number))
        })
        .collect()
    }
  }
}
