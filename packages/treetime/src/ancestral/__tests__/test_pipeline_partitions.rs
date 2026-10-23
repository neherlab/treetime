#[cfg(test)]
mod tests {
  use crate::ancestral::params::MethodAncestral;
  use eyre::Report;
  use helpers::{assert_same_alignment, assert_same_partition, column_blocks, reconstruct, total_log_lh};
  use itertools::{Itertools, izip};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::slice;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::fitch( MethodAncestral::Parsimony, None       )]
  #[case::dense( MethodAncestral::Marginal,  Some(true) )]
  #[case::sparse(MethodAncestral::Marginal,  Some(false))]
  #[trace]
  fn test_ancestral_two_alignments_equal_two_separate_runs(
    #[case] method: MethodAncestral,
    #[case] dense: Option<bool>,
  ) -> Result<(), Report> {
    let blocks = column_blocks()?;
    let (joint, joint_tracks) = reconstruct(method, dense, &blocks)?;
    let separate = blocks
      .iter()
      .map(|block| reconstruct(method, dense, slice::from_ref(block)))
      .collect::<Result<Vec<_>, Report>>()?;

    assert_eq!(
      (vec![1, 2], vec![(1, 1), (2, 2)]),
      (
        joint.alignments.iter().map(|alignment| alignment.number).collect_vec(),
        joint
          .partitions
          .iter()
          .map(|partition| (partition.number, partition.alignment))
          .collect_vec()
      )
    );
    for (joint_alignment, joint_partition, (single, _)) in izip!(&joint.alignments, &joint.partitions, &separate) {
      let ([single_alignment], [single_partition]) = (single.alignments.as_slice(), single.partitions.as_slice()) else {
        panic!("one alignment gives one alignment result and one partition result");
      };
      assert_same_alignment(joint_alignment, single_alignment);
      assert_same_partition(joint_partition, single_partition);
    }

    let node_count = joint.graph.num_nodes();
    let expected_tracks = [vec![1; node_count], vec![2; node_count]].concat();
    assert_eq!(expected_tracks, joint_tracks);

    let separate_total = separate
      .iter()
      .map(|(single, _)| total_log_lh(&single.partitions))
      .sum::<Option<f64>>();
    match (total_log_lh(&joint.partitions), separate_total) {
      (Some(joint_total), Some(separate_total)) => {
        pretty_assert_abs_diff_eq!(separate_total, joint_total, epsilon = 1e-10);
      },
      (joint_total, separate_total) => {
        assert_eq!(
          (None, None, MethodAncestral::Parsimony),
          (joint_total, separate_total, method)
        );
      },
    }
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::ancestral::pipeline::{AncestralOutput, run};
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModel;
    use crate::partition::marginal::sample::SampleMode;
    use crate::partition::result::{AlignmentResult, PartitionResult};
    use crate::progress::NoopProgress;
    use crate::seq::sink::SeqStream;
    use crate::seq::sink::{SeqItem, SeqSink, SeqTrack};
    use crate::test_utils::{TestDataset, ancestral_input_for_tree, complete_sequence_source, get_common_length};
    use eyre::Report;
    use pretty_assertions::assert_eq;
    use std::ops::Range;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, LogLh, Seq};
    use treetime_utils::pretty_assert_abs_diff_eq;

    pub(super) fn column_blocks() -> Result<Vec<Block>, Report> {
      let records = TestDataset::Flu20.alignment()?;
      let length = get_common_length(&records)?;
      let split = length.div_ceil(2);
      let block = |path: &str, columns: Range<usize>| Block {
        path: path.to_owned(),
        records: records
          .iter()
          .map(|record| AlignmentRecord {
            name: record.name.clone(),
            seq: Seq::from(&record.seq.as_slice()[columns.clone()]),
          })
          .collect(),
      };
      Ok(vec![
        block("first.fasta", 0..split),
        block("second.fasta", split..length),
      ])
    }

    #[derive(Clone)]
    pub(super) struct Block {
      path: String,
      records: Vec<AlignmentRecord>,
    }

    pub(super) fn reconstruct(
      method: MethodAncestral,
      dense: Option<bool>,
      blocks: &[Block],
    ) -> Result<(AncestralOutput, Vec<usize>), Report> {
      let alphabet = Alphabet::default();
      let parse = nwk_read(TestDataset::Flu20.tree().as_bytes())?;
      let names = parse.names.clone();
      let sources = blocks
        .iter()
        .map(|block| complete_sequence_source(&parse.graph, &names, block.records.clone(), &alphabet, &block.path))
        .collect::<Result<Vec<_>, Report>>()?;
      let params = AncestralParams {
        method,
        dense,
        include_leaves: true,
        report_ambiguous: true,
        impute_missing_data: false,
        seed: 0,
        sample_from_profile: SampleMode::Argmax,
      };
      let mut sink = TrackSink::default();
      let output = run(
        &params,
        ancestral_input_for_tree(parse, sources, &GtrModel::Infer)?,
        Some(SeqStream::new(&mut sink, &names)),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      Ok((output, sink.tracks))
    }

    pub(super) fn assert_same_alignment(joint: &AlignmentResult, single: &AlignmentResult) {
      assert_eq!(
        (&single.path, single.length, &single.mask, single.indel_rate),
        (&joint.path, joint.length, &joint.mask, joint.indel_rate)
      );
      let joint_sequences = joint.sequences.as_ref().expect("ancestral returns sequences");
      let single_sequences = single.sequences.as_ref().expect("ancestral returns sequences");
      assert_eq!(single_sequences.root_sequence, joint_sequences.root_sequence);
      assert_eq!(single_sequences.edge_mutations, joint_sequences.edge_mutations);
      assert_eq!(
        single_sequences.edge_mutation_counts,
        joint_sequences.edge_mutation_counts
      );
    }

    pub(super) fn assert_same_partition(joint: &PartitionResult, single: &PartitionResult) {
      assert_eq!(
        (single.sites, single.model_name, &single.name),
        (joint.sites, joint.model_name, &joint.name)
      );
      let model = |result: &PartitionResult| {
        result
          .gtr
          .as_ref()
          .map(|gtr| (gtr.mu(), gtr.W().clone(), gtr.pi().clone()))
      };
      assert_eq!(model(single), model(joint));
      match (single.log_lh, joint.log_lh) {
        (Some(single), Some(joint)) => pretty_assert_abs_diff_eq!(single.value(), joint.value(), epsilon = 1e-10),
        (single, joint) => assert_eq!((None, None), (single, joint)),
      }
    }

    pub(super) fn total_log_lh(results: &[PartitionResult]) -> Option<f64> {
      results.iter().map(|result| result.log_lh.map(LogLh::value)).sum()
    }

    #[derive(Default)]
    struct TrackSink {
      tracks: Vec<usize>,
    }

    impl SeqSink for TrackSink {
      fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
        let SeqTrack::Nuc(number) = item.track else {
          panic!("ancestral reconstruction streams only nucleotide tracks");
        };
        self.tracks.push(number);
        Ok(())
      }
    }
  }
}
