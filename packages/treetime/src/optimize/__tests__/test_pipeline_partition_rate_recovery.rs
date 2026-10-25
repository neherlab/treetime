#[cfg(test)]
mod tests {
  use crate::test_utils::{recovery_check, require_passed, run_replicates};
  use eyre::Report;
  use helpers::{BASE_SEED, FILE_COUNT, MARGIN, REPLICATE_COUNT, SIMULATED_RATES, replicate_rates};
  use itertools::Itertools;
  use ndarray::{Array1, Array2};
  use treetime_simulation::FALSE_FAIL_RATE;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::make_internal_report;

  #[test]
  fn test_sim_optimize_partition_rates_recover_the_simulated_relative_rates() -> Result<(), Report> {
    let simulated = Array1::from_iter(SIMULATED_RATES);
    let targets = &simulated
      / simulated
        .mean()
        .ok_or_else(|| make_internal_report!("No simulated rates"))?;

    let estimates = run_replicates(BASE_SEED, REPLICATE_COUNT, replicate_rates)?;

    let estimates = Array2::from_shape_vec((REPLICATE_COUNT, FILE_COUNT), estimates.concat())?;
    let log_ratios = (&estimates / &targets).ln();
    let checks = log_ratios
      .columns()
      .into_iter()
      .zip(SIMULATED_RATES)
      .map(|(errors, rate)| {
        recovery_check(
          &format!("log ratio of the estimated and the simulated relative rate {rate}"),
          &errors.to_vec(),
          MARGIN,
          FALSE_FAIL_RATE / FILE_COUNT.to_f64_exact(),
        )
      })
      .try_collect::<_, Vec<_>, _>()?;
    require_passed(&checks)
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
    use crate::partition::split::SequenceInput;
    use crate::progress::NoopProgress;
    use crate::seq::alignment::SequenceSource;
    use crate::test_utils::{TestDataset, simulated_sequence_source};
    use eyre::Report;
    use itertools::Itertools;
    use treetime_graph::key_map::EdgeMap;

    use std::path::PathBuf;

    use treetime_io::nwk::nwk_read;
    use treetime_simulation::evolve::evolve_sequences;
    use treetime_simulation::substitution_model::SubstitutionModel;
    use treetime_utils::sync::random::get_random_number_generator;
    use treetime_utils::{make_internal_report, make_report};

    pub(super) const SIMULATED_RATES: [f64; 3] = [0.3, 1.0, 2.0];
    pub(super) const FILE_COUNT: usize = SIMULATED_RATES.len();
    pub(super) const REPLICATE_COUNT: usize = 10275;
    pub(super) const BASE_SEED: u64 = 1;
    pub(super) const MARGIN: f64 = 0.03;
    const SITES_PER_FILE: usize = 3000;

    pub(super) fn replicate_rates(seed: u64) -> Result<Vec<f64>, Report> {
      let parse = nwk_read(TestDataset::Zika86.tree().as_bytes())?;
      let durations: EdgeMap<f64> = parse
        .branch_lengths
        .iter()
        .map(|(edge, &length)| {
          length
            .map(|length| (edge, length))
            .ok_or_else(|| make_report!("Edge {edge} of the Zika tree has no branch length"))
        })
        .try_collect()?;
      let mut rng = get_random_number_generator(seed);
      let sources = SIMULATED_RATES
        .iter()
        .enumerate()
        .map(|(file, &rate)| {
          let model = SubstitutionModel::jc69(rate)?;
          let sequences = evolve_sequences(&parse.graph, &durations, &model, SITES_PER_FILE, &mut rng)?;
          let source = simulated_sequence_source(&parse.graph, &parse.names, &sequences.sequences)?;
          Ok(SequenceSource {
            path: PathBuf::from(format!("{}.fasta", file + 1)),
            ..source
          })
        })
        .collect::<Result<Vec<_>, Report>>()?;
      let params = OptimizeParams {
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
        names: parse.names,
        sequences: SequenceInput::whole_alignments(sources, &GtrModel::JC69(JC69Params::default()))?,
        branch_lengths: parse.branch_lengths,
      };
      let output = run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress)?;
      output
        .partitions
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
