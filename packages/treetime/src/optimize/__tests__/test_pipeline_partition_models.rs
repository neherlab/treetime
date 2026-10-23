#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{GtrModel, JC69Params};
  use eyre::Report;
  use helpers::{block_rates, optimize};
  use pretty_assertions::assert_eq;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_optimize_normalizes_an_inferred_and_a_named_partition_before_the_first_iteration() -> Result<(), Report> {
    let models = [GtrModel::Infer, GtrModel::JC69(JC69Params::default())];

    let rates = block_rates(&optimize(&models, 0)?);

    let weighted_mean = rates.iter().map(|(sites, mu)| sites.to_f64_exact() * mu).sum::<f64>()
      / rates.iter().map(|(sites, _)| sites.to_f64_exact()).sum::<f64>();
    pretty_assert_ulps_eq!(1.0, weighted_mean, max_ulps = 4);
    assert!(
      (rates[1].1 - 1.0).abs() > 1e-6,
      "the named partition shares the normalization, got {rates:?}"
    );
    Ok(())
  }

  #[test]
  fn test_optimize_without_an_inferred_partition_keeps_the_rates_of_the_named_models() -> Result<(), Report> {
    let models = [
      GtrModel::JC69(JC69Params::default()),
      GtrModel::JC69(JC69Params::default()),
    ];

    let rates = block_rates(&optimize(&models, 0)?);

    assert_eq!(vec![1.0, 1.0], rates.iter().map(|(_, mu)| *mu).collect::<Vec<_>>());
    Ok(())
  }

  #[test]
  fn test_optimize_one_partition_keeps_the_rate_of_its_named_model() -> Result<(), Report> {
    let models = [GtrModel::JC69(JC69Params {
      mu: 2.0,
      ..JC69Params::default()
    })];

    let rates = block_rates(&optimize(&models, 1)?);

    assert_eq!(vec![2.0], rates.iter().map(|(_, mu)| *mu).collect::<Vec<_>>());
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModel;
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::pipeline::{OptimizeInput, OptimizeOutput, OptimizeParams, run};
    use crate::partition::resolve::{PartitionEntry, resolve_partitions};
    use crate::partition::split::SequenceInput;
    use crate::progress::NoopProgress;
    use crate::seq::alignment::SequenceSource;
    use crate::test_utils::{Columns, TestDataset, leaf_seq_inputs, sequence_source, split_at_gap_free_column};
    use eyre::Report;
    use itertools::{Itertools, izip};
    use std::path::PathBuf;
    use treetime_io::nwk::nwk_read;

    pub(super) fn optimize(models: &[GtrModel], max_iter: usize) -> Result<OptimizeOutput, Report> {
      let dataset = TestDataset::Flu20;
      let parse = nwk_read(dataset.tree().as_bytes())?;
      let records = dataset.alignment_columns(Columns::Leading(400))?;
      let blocks = if models.len() == 1 {
        vec![records]
      } else {
        split_at_gap_free_column(&records)?.into_iter().collect_vec()
      };
      let alignments = izip!(1.., blocks)
        .map(|(number, records)| {
          let source = sequence_source(
            Alphabet::default(),
            leaf_seq_inputs(&parse.graph, &parse.names, records),
          )?;
          Ok(SequenceSource {
            path: PathBuf::from(format!("{number}.fasta")),
            ..source
          })
        })
        .collect::<Result<Vec<_>, Report>>()?;
      let entries = izip!(1.., models)
        .map(|(alignment, model)| PartitionEntry {
          name: None,
          alignment: Some(alignment),
          sites: None,
          model: model.clone(),
        })
        .collect_vec();
      let lengths = alignments.iter().map(|alignment| alignment.length).collect_vec();
      let partitions = resolve_partitions(&entries, &lengths, &GtrModel::Infer, None)?.plans;
      let params = OptimizeParams {
        dense: Some(true),
        initial_guess: InitialGuessMode::Never,
        reroot_spec: None,
        loop_settings: OptimizeLoopSettings {
          max_iter,
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
        sequences: SequenceInput { alignments, partitions },
        branch_lengths: parse.branch_lengths,
      };
      Ok(run(&params, input, &NoopCancel, &NoopProgress, &NoopProgress)?)
    }

    pub(super) fn block_rates(output: &OptimizeOutput) -> Vec<(usize, f64)> {
      output
        .partitions
        .iter()
        .map(|partition| {
          let gtr = partition.gtr.as_ref().expect("optimize fits a model");
          (partition.sites, gtr.mu())
        })
        .collect()
    }
  }
}
