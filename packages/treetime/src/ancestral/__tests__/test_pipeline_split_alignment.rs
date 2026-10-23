#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{GtrModel, JC69Params};
  use crate::partition::create::Representation;
  use crate::test_utils::{Columns, TestDataset};
  use eyre::Report;
  use helpers::{MixedModels, Reconstructed, codon_entries, columns, reconstruct};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::zika_dense(  TestDataset::Zika20,  Representation::Dense )]
  #[case::zika_sparse( TestDataset::Zika20,  Representation::Sparse)]
  #[case::flu_dense(   TestDataset::Flu20,   Representation::Dense )]
  #[case::flu_sparse(  TestDataset::Flu20,   Representation::Sparse)]
  #[case::ebola_dense( TestDataset::Ebola20, Representation::Dense )]
  #[case::ebola_sparse(TestDataset::Ebola20, Representation::Sparse)]
  #[trace]
  fn test_ancestral_codon_positions_give_the_sequences_and_mutations_of_one_alignment(
    #[case] dataset: TestDataset,
    #[case] representation: Representation,
  ) -> Result<(), Report> {
    let records = dataset.alignment_columns(Columns::WholeCodons)?;
    let jc69 = GtrModel::JC69(JC69Params::default());
    let whole = reconstruct(dataset, &records, representation, &[])?;
    let split = reconstruct(dataset, &records, representation, &codon_entries(&records, &jc69)?)?;

    let Reconstructed {
      emitted,
      root_sequence,
      edge_mutations,
      edge_mutation_counts,
      log_lh,
      ..
    } = whole;
    assert_eq!(
      (emitted, root_sequence, edge_mutations, edge_mutation_counts),
      (split.emitted, split.root_sequence, split.edge_mutations, split.edge_mutation_counts)
    );
    pretty_assert_abs_diff_eq!(log_lh.iter().sum::<f64>(), split.log_lh.iter().sum::<f64>(), epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense_first_codon_positions(  Representation::Dense,  0)]
  #[case::dense_other_codon_positions(  Representation::Dense,  1)]
  #[case::sparse_first_codon_positions( Representation::Sparse, 0)]
  #[case::sparse_other_codon_positions( Representation::Sparse, 1)]
  #[trace]
  fn test_ancestral_split_alignment_with_two_models_gives_each_partition_the_result_of_its_separate_run(
    #[case] representation: Representation,
    #[case] partition: usize,
  ) -> Result<(), Report> {
    let case = MixedModels::new()?;
    let split = case.split(representation)?;
    let single = case.separate(partition, representation)?;

    assert_eq!(
      (&single.models[0], single.emitted),
      (&split.models[partition], columns(&split.emitted, &case.columns[partition]))
    );
    pretty_assert_abs_diff_eq!(single.log_lh[0], split.log_lh[partition], epsilon = 1e-8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( Representation::Dense )]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_ancestral_split_alignment_with_two_models_sums_the_log_likelihoods_of_the_separate_runs(
    #[case] representation: Representation,
  ) -> Result<(), Report> {
    let case = MixedModels::new()?;
    let split = case.split(representation)?;
    let separate = [case.separate(0, representation)?, case.separate(1, representation)?];

    pretty_assert_abs_diff_eq!(
      separate.iter().map(|single| single.log_lh[0]).sum::<f64>(),
      split.log_lh.iter().sum::<f64>(),
      epsilon = 1e-8
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::ancestral::pipeline::run;
    use crate::ancestral::plan::{ReconstructionOptions, reconstruct_marginal};
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::cancel::NoopCancel;
    use crate::error::OperationError;
    use crate::gtr::get_gtr::{GtrModel, HKY85Params, JC69Params};
    use crate::gtr::gtr::GTR;
    use crate::partition::alignment_sequences::PartitionSequences;
    use crate::partition::create::{Representation, build_marginal_partition};
    use crate::partition::inferred_rate::InferredRate;
    use crate::partition::marginal::posterior_states::AmbiguityThreshold;
    use crate::partition::marginal::sample::SampleMode;
    use crate::partition::resolve::{PartitionEntry, resolve_partitions};
    use crate::partition::sites::SiteRange;
    use crate::partition::split::SequenceInput;
    use crate::progress::NoopProgress;
    use crate::seq::alignment::{AncestralInput, EdgeSeqInput};
    use crate::seq::mutation::Mutation;
    use crate::seq::sink::{AmbiguousSites, SeqItem, SeqSink, SeqStream, SeqTrack};
    use crate::test_utils::{Columns, TestDataset, complete_sequence_source, get_common_length};
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::{Array1, Array2};
    use std::collections::BTreeMap;
    use std::str::FromStr;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, LogLh, Seq};
    use treetime_utils::sync::random::get_random_number_generator;

    const AMBIGUITY_THRESHOLD: f64 = 0.2;

    pub(super) struct Reconstructed {
      pub(super) emitted: BTreeMap<GraphNodeKey, EmittedSequence>,
      pub(super) root_sequence: Seq,
      pub(super) edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
      pub(super) edge_mutation_counts: BTreeMap<GraphEdgeKey, usize>,
      pub(super) models: Vec<Model>,
      pub(super) log_lh: Vec<f64>,
    }

    #[derive(Clone, Debug, PartialEq)]
    pub(super) struct EmittedSequence {
      pub(super) track: usize,
      pub(super) emitted: bool,
      pub(super) seq: Seq,
      pub(super) ambiguous_sites: AmbiguousSites,
    }

    #[derive(Clone, Debug, PartialEq)]
    pub(super) struct Model {
      w: Array2<f64>,
      pi: Array1<f64>,
      mu: f64,
    }

    pub(super) struct MixedModels {
      dataset: TestDataset,
      records: Vec<AlignmentRecord>,
      models: [GtrModel; 2],
      pub(super) columns: [Vec<usize>; 2],
    }

    impl MixedModels {
      pub(super) fn new() -> Result<Self, Report> {
        let dataset = TestDataset::Flu20;
        let hky85 = GtrModel::HKY85(HKY85Params {
          kappa: 4.0,
          ..HKY85Params::default()
        });
        Ok(Self {
          dataset,
          records: dataset.alignment_columns(Columns::Leading(900))?,
          models: [GtrModel::JC69(JC69Params::default()), hky85],
          columns: [
            (0..900).step_by(3).collect(),
            (0..900).filter(|column| column % 3 != 0).collect(),
          ],
        })
      }

      pub(super) fn split(&self, representation: Representation) -> Result<Reconstructed, Report> {
        let [jc69, hky85] = &self.models;
        let entries = [sites_entry("1-900\\3", jc69)?, sites_entry("2-900\\3,3-900\\3", hky85)?];
        reconstruct(self.dataset, &self.records, representation, &entries)
      }

      pub(super) fn separate(&self, partition: usize, representation: Representation) -> Result<Reconstructed, Report> {
        separate_run(
          self.dataset,
          &self.records,
          &self.columns[partition],
          &self.models[partition],
          representation,
        )
      }
    }

    pub(super) fn codon_entries(records: &[AlignmentRecord], model: &GtrModel) -> Result<Vec<PartitionEntry>, Report> {
      let length = get_common_length(records)?;
      (1..=3)
        .map(|position| sites_entry(&format!("{position}-{length}\\3"), model))
        .collect()
    }

    fn sites_entry(sites: &str, model: &GtrModel) -> Result<PartitionEntry, Report> {
      Ok(PartitionEntry {
        name: None,
        alignment: Some(1),
        sites: Some(sites.split(',').map(SiteRange::from_str).try_collect()?),
        model: model.clone(),
      })
    }

    pub(super) fn reconstruct(
      dataset: TestDataset,
      records: &[AlignmentRecord],
      representation: Representation,
      entries: &[PartitionEntry],
    ) -> Result<Reconstructed, Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read(dataset.tree().as_bytes())?;
      let names = parse.names.clone();
      let source = complete_sequence_source(&parse.graph, &names, records.to_vec(), &alphabet, "aln.fasta")?;
      let default_model = entries
        .first()
        .map_or_else(|| GtrModel::JC69(JC69Params::default()), |entry| entry.model.clone());
      let partitions = resolve_partitions(entries, &[source.length], &default_model, Some(alphabet.name()))?.plans;
      let input = AncestralInput {
        edges: parse
          .branch_lengths
          .into_iter()
          .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
          .collect(),
        graph: parse.graph,
        sequences: SequenceInput {
          alignments: vec![source],
          partitions,
        },
      };
      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        dense: Some(matches!(representation, Representation::Dense)),
        include_leaves: true,
        report_ambiguous: true,
        impute_missing_data: false,
        seed: 0,
        sample_from_profile: SampleMode::Argmax,
      };
      let mut sink = AmbiguitySink::new()?;
      let output = run(
        &params,
        input,
        Some(SeqStream::new(&mut sink, &names)),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let [alignment] = output.alignments.as_slice() else {
        panic!("one alignment gives one alignment result");
      };
      let sequences = alignment.sequences.clone().expect("ancestral returns sequences");
      Ok(Reconstructed {
        emitted: sink.emitted,
        root_sequence: sequences.root_sequence,
        edge_mutations: sequences.edge_mutations,
        edge_mutation_counts: sequences.edge_mutation_counts,
        models: output
          .partitions
          .iter()
          .map(|partition| model_of(partition.gtr.as_ref().expect("a marginal run fits a model")))
          .collect(),
        log_lh: output
          .partitions
          .iter()
          .map(|partition| {
            partition
              .log_lh
              .map(LogLh::value)
              .expect("a marginal run has a likelihood")
          })
          .collect(),
      })
    }

    fn separate_run(
      dataset: TestDataset,
      records: &[AlignmentRecord],
      columns_of: &[usize],
      model: &GtrModel,
      representation: Representation,
    ) -> Result<Reconstructed, Report> {
      let alphabet = Alphabet::default();
      let alignment_length = get_common_length(records)?;
      let parse = nwk_read(dataset.tree().as_bytes())?;
      let subset = records
        .iter()
        .map(|record| AlignmentRecord {
          name: record.name.clone(),
          seq: columns_of.iter().map(|&column| record.seq[column]).collect(),
        })
        .collect_vec();
      let source = complete_sequence_source(&parse.graph, &parse.names, subset, &alphabet, "subset.fasta")?;
      let branch_lengths = branch_lengths_or_zero(&parse.branch_lengths);
      let reconstruction = build_marginal_partition(
        representation,
        model,
        InferredRate::Unit,
        &parse.graph,
        source.alphabet,
        source.nodes,
        &branch_lengths,
        &NoopProgress,
      )?
      .with_total_length(alignment_length);
      let partition = reconstruct_marginal(
        &parse.graph,
        reconstruction,
        &branch_lengths,
        &ReconstructionOptions::new(false, SampleMode::Argmax),
        &mut get_random_number_generator(0),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      let mut sink = AmbiguitySink::new()?;
      let mutations = partition
        .stream_sequences(
          &parse.graph,
          SeqTrack::Nuc(1),
          true,
          true,
          Some(SeqStream::new(&mut sink, &parse.names)),
        )
        .map_err(OperationError::into_report)?;
      let log_lh = partition
        .log_lh(&parse.graph)?
        .expect("a marginal run has a likelihood");
      Ok(Reconstructed {
        emitted: sink.emitted,
        root_sequence: mutations.root_sequence,
        edge_mutations: mutations.edge_mutations,
        edge_mutation_counts: BTreeMap::new(),
        models: vec![model_of(partition.gtr().expect("a marginal run fits a model"))],
        log_lh: vec![log_lh.value()],
      })
    }

    pub(super) fn columns(
      emitted: &BTreeMap<GraphNodeKey, EmittedSequence>,
      columns_of: &[usize],
    ) -> BTreeMap<GraphNodeKey, EmittedSequence> {
      emitted
        .iter()
        .map(|(&key, sequence)| {
          let ambiguous_sites = columns_of
            .iter()
            .enumerate()
            .filter_map(|(local, column)| sequence.ambiguous_sites.get(column).map(|&code| (local, code)))
            .collect();
          let projected = EmittedSequence {
            track: sequence.track,
            emitted: sequence.emitted,
            seq: columns_of.iter().map(|&column| sequence.seq[column]).collect(),
            ambiguous_sites,
          };
          (key, projected)
        })
        .collect()
    }

    fn model_of(gtr: &GTR) -> Model {
      Model {
        w: gtr.W().clone(),
        pi: gtr.pi().clone(),
        mu: gtr.mu(),
      }
    }

    struct AmbiguitySink {
      threshold: AmbiguityThreshold,
      emitted: BTreeMap<GraphNodeKey, EmittedSequence>,
    }

    impl AmbiguitySink {
      fn new() -> Result<Self, Report> {
        Ok(Self {
          threshold: AmbiguityThreshold::new(AMBIGUITY_THRESHOLD)?,
          emitted: BTreeMap::new(),
        })
      }
    }

    impl SeqSink for AmbiguitySink {
      fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
        let SeqTrack::Nuc(track) = item.track else {
          panic!("ancestral reconstruction streams only nucleotide tracks");
        };
        let sequence = EmittedSequence {
          track,
          emitted: item.emitted,
          seq: item.seq.clone(),
          ambiguous_sites: item.ambiguous_sites.clone(),
        };
        assert!(
          self.emitted.insert(item.key, sequence).is_none(),
          "each node is emitted once per alignment"
        );
        Ok(())
      }

      fn ambiguity(&self) -> Option<AmbiguityThreshold> {
        Some(self.threshold)
      }
    }
  }
}
