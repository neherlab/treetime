#[cfg(test)]
mod tests {
  use crate::ancestral::plan::ReconstructionMethod;
  use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
  use crate::gtr::get_gtr::{GtrModel, GtrModelName, HKY85Params, JC69Params};
  use crate::partition::create::Representation;
  use crate::test_utils::RecordingLog;
  use eyre::Report;
  use helpers::{
    Fixture, JC69TwoTaxaSite, TWO_TAXA_LENGTH, codon_entry, entry, jc69_two_taxa_site_likelihoods, log_lh_by_partition,
    rates_by_partition, reconstruct, run_pipeline, sites_entry,
  };
  use itertools::Itertools;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( Representation::Dense )]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_ancestral_partitions_of_a_split_alignment_floor_zero_length_branches_at_the_alignment_length(
    #[case] representation: Representation,
  ) -> Result<(), Report> {
    let fixture = Fixture::two_taxa()?;
    let jc69 = GtrModel::JC69(JC69Params::default());
    let entries = vec![sites_entry("1-2", jc69.clone()), sites_entry("3-4", jc69)];

    let actual = log_lh_by_partition(&fixture, &entries, ReconstructionMethod::Marginal(representation))?;

    let floor = MIN_BRANCH_LENGTH_FRACTION / TWO_TAXA_LENGTH.to_f64_exact();
    let JC69TwoTaxaSite { same, different } = jc69_two_taxa_site_likelihoods(2.0 * floor);
    let expected = array![2.0 * same.ln(), same.ln() + different.ln()];
    pretty_assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( Representation::Dense )]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_ancestral_inferred_rates_of_a_split_alignment_have_a_site_weighted_mean_of_one(
    #[case] representation: Representation,
  ) -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let mut entries = (1..=3).map(|position| codon_entry(position, GtrModel::Infer)).collect_vec();
    entries.push(entry(2, None, GtrModel::Infer));

    let partitions = reconstruct(&fixture, &entries, ReconstructionMethod::Marginal(representation))?;

    let rates = rates_by_partition(&partitions);
    let split = &rates[..3];
    let weighted_mean = split.iter().map(|(sites, mu)| sites.to_f64_exact() * mu).sum::<f64>()
      / split.iter().map(|(sites, _)| sites.to_f64_exact()).sum::<f64>();
    pretty_assert_ulps_eq!(1.0, weighted_mean, max_ulps = 4);
    assert!(
      split.iter().map(|(_, mu)| mu.to_bits()).all_unique(),
      "the codon positions fit different rates, got {split:?}"
    );
    pretty_assert_ulps_eq!(1.0, rates[3].1, max_ulps = 0);
    Ok(())
  }

  #[test]
  fn test_ancestral_named_partitions_of_a_split_alignment_keep_the_rate_of_their_model() -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let entries = vec![
      codon_entry(1, GtrModel::Infer),
      codon_entry(2, GtrModel::Infer),
      codon_entry(
        3,
        GtrModel::JC69(JC69Params {
          mu: 2.0,
          ..JC69Params::default()
        }),
      ),
      entry(2, None, GtrModel::Infer),
    ];

    let partitions = reconstruct(
      &fixture,
      &entries,
      ReconstructionMethod::Marginal(Representation::Dense),
    )?;

    let rates = rates_by_partition(&partitions);
    let inferred = &rates[..2];
    let weighted_mean = inferred
      .iter()
      .map(|(sites, mu)| sites.to_f64_exact() * mu)
      .sum::<f64>()
      / inferred.iter().map(|(sites, _)| sites.to_f64_exact()).sum::<f64>();
    pretty_assert_ulps_eq!(1.0, weighted_mean, max_ulps = 4);
    assert_eq!((2.0, 1.0), (rates[2].1, rates[3].1));
    Ok(())
  }

  #[test]
  fn test_ancestral_partition_results_carry_the_model_of_their_entry() -> Result<(), Report> {
    let fixture = Fixture::new()?;
    let hky85 = GtrModel::HKY85(HKY85Params {
      kappa: 4.0,
      ..HKY85Params::default()
    });
    let entries = vec![entry(1, None, GtrModel::Infer), entry(2, None, hky85.clone())];

    let output = run_pipeline(&fixture, &entries, &RecordingLog::default())?;

    let summary = output
      .partitions
      .iter()
      .map(|partition| (partition.number, partition.alignment, partition.model_name))
      .collect_vec();
    assert_eq!(vec![(1, 1, GtrModelName::Infer), (2, 2, GtrModelName::HKY85)], summary);
    let expected = hky85.build(fixture.alphabet.name())?;
    let actual = output.partitions[1].gtr.as_ref().expect("a marginal run fits a model");
    assert_eq!(
      (expected.W().clone(), expected.pi().clone(), expected.mu()),
      (actual.W().clone(), actual.pi().clone(), actual.mu())
    );
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::ancestral::pipeline::{AncestralOutput, run};
    use crate::ancestral::plan::{
      ReconstructedPartition, ReconstructionMethod, ReconstructionOptions, reconstruct_partitions,
    };
    use crate::branch_lengths::branch_lengths_or_zero;
    use crate::cancel::NoopCancel;
    use crate::error::OperationError;
    use crate::gtr::get_gtr::GtrModel;
    use crate::partition::marginal::sample::SampleMode;
    use crate::partition::resolve::{PartitionEntry, resolve_partitions};
    use crate::partition::sites::SiteRange;
    use crate::partition::split::SequenceInput;
    use crate::progress::{LogSink, NoopProgress};
    use crate::seq::alignment::{AncestralInput, EdgeSeqInput};
    use crate::test_utils::{TestDataset, complete_sequence_source, get_common_length};
    use eyre::Report;
    use itertools::Itertools;
    use ndarray::Array1;
    use std::str::FromStr;
    use treetime_graph::graph::Graph;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AlignmentRecord, Seq};
    use treetime_utils::sync::random::get_random_number_generator;

    pub(super) const CODON_BLOCK: usize = 900;

    pub(super) const TWO_TAXA_LENGTH: usize = 4;

    pub(super) struct Fixture {
      pub(super) alphabet: Alphabet,
      tree: &'static str,
      alignments: Vec<(&'static str, Vec<AlignmentRecord>)>,
    }

    impl Fixture {
      pub(super) fn new() -> Result<Self, Report> {
        let records = TestDataset::Flu20.alignment()?;
        let length = get_common_length(&records)?;
        let block = |first: usize, last: usize| {
          records
            .iter()
            .map(|record| AlignmentRecord {
              name: record.name.clone(),
              seq: Seq::from(&record.seq.as_slice()[first..last]),
            })
            .collect_vec()
        };
        Ok(Self {
          alphabet: Alphabet::default(),
          tree: TestDataset::Flu20.tree(),
          alignments: vec![
            ("first.fasta", block(0, CODON_BLOCK)),
            ("second.fasta", block(CODON_BLOCK, length)),
          ],
        })
      }

      pub(super) fn two_taxa() -> Result<Self, Report> {
        let record = |name: &str, seq: &str| -> Result<AlignmentRecord, Report> {
          Ok(AlignmentRecord {
            name: name.to_owned(),
            seq: Seq::try_from_str(seq)?,
          })
        };
        Ok(Self {
          alphabet: Alphabet::default(),
          tree: "(A:0,B:0)root;",
          alignments: vec![("two_taxa.fasta", vec![record("A", "ACGT")?, record("B", "ACGA")?])],
        })
      }

      fn input(&self, entries: &[PartitionEntry]) -> Result<AncestralInput, Report> {
        let parse = nwk_read(self.tree.as_bytes())?;
        let alignments = self
          .alignments
          .iter()
          .map(|(path, records)| {
            complete_sequence_source(&parse.graph, &parse.names, records.clone(), &self.alphabet, path)
          })
          .collect::<Result<Vec<_>, Report>>()?;
        let lengths = alignments.iter().map(|alignment| alignment.length).collect_vec();
        let partitions = resolve_partitions(entries, &lengths, &GtrModel::Infer, Some(self.alphabet.name()))?.plans;
        Ok(AncestralInput {
          edges: parse
            .branch_lengths
            .into_iter()
            .map(|(key, branch_length)| (key, EdgeSeqInput { branch_length }))
            .collect(),
          graph: parse.graph,
          sequences: SequenceInput { alignments, partitions },
        })
      }
    }

    pub(super) fn codon_entry(position: usize, model: GtrModel) -> PartitionEntry {
      let range = SiteRange::from_str(&format!("{position}-{CODON_BLOCK}\\3")).expect("the codon range is valid");
      entry(1, Some(vec![range]), model)
    }

    pub(super) fn sites_entry(sites: &str, model: GtrModel) -> PartitionEntry {
      let range = SiteRange::from_str(sites).expect("the site range is valid");
      entry(1, Some(vec![range]), model)
    }

    pub(super) fn entry(alignment: usize, sites: Option<Vec<SiteRange>>, model: GtrModel) -> PartitionEntry {
      PartitionEntry {
        name: None,
        alignment: Some(alignment),
        sites,
        model,
      }
    }

    pub(super) fn reconstruct(
      fixture: &Fixture,
      entries: &[PartitionEntry],
      method: ReconstructionMethod,
    ) -> Result<Vec<ReconstructedPartition>, Report> {
      Ok(reconstruct_with_graph(fixture, entries, method)?.1)
    }

    pub(super) fn log_lh_by_partition(
      fixture: &Fixture,
      entries: &[PartitionEntry],
      method: ReconstructionMethod,
    ) -> Result<Array1<f64>, Report> {
      let (graph, partitions) = reconstruct_with_graph(fixture, entries, method)?;
      partitions
        .iter()
        .map(|partition| {
          let log_lh = partition
            .partition
            .log_lh(&graph)?
            .expect("a marginal partition has a likelihood");
          Ok(log_lh.value())
        })
        .collect()
    }

    pub(super) struct JC69TwoTaxaSite {
      pub(super) same: f64,
      pub(super) different: f64,
    }

    pub(super) fn jc69_two_taxa_site_likelihoods(path_length: f64) -> JC69TwoTaxaSite {
      let decay = (-4.0 / 3.0 * path_length).exp();
      let change = -(-4.0 / 3.0 * path_length).exp_m1();
      JC69TwoTaxaSite {
        same: 0.25 * (0.25 + 0.75 * decay),
        different: 0.25 * 0.25 * change,
      }
    }

    fn reconstruct_with_graph(
      fixture: &Fixture,
      entries: &[PartitionEntry],
      method: ReconstructionMethod,
    ) -> Result<(Graph, Vec<ReconstructedPartition>), Report> {
      let input = fixture.input(entries)?;
      let branch_lengths = branch_lengths_or_zero(&input.branch_lengths());
      let AncestralInput { graph, sequences, .. } = input;
      let partitions = reconstruct_partitions(
        &graph,
        method,
        sequences.split()?,
        &branch_lengths,
        &ReconstructionOptions::new(false, SampleMode::Argmax),
        &mut get_random_number_generator(0),
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?;
      Ok((graph, partitions))
    }

    pub(super) fn run_pipeline(
      fixture: &Fixture,
      entries: &[PartitionEntry],
      log: &dyn LogSink,
    ) -> Result<AncestralOutput, OperationError> {
      let input = fixture.input(entries).map_err(OperationError::InvalidInput)?;
      let params = AncestralParams {
        method: MethodAncestral::Marginal,
        dense: Some(true),
        include_leaves: true,
        report_ambiguous: true,
        impute_missing_data: false,
        seed: 0,
        sample_from_profile: SampleMode::Argmax,
      };
      run(&params, input, None, &NoopCancel, &NoopProgress, log)
    }

    pub(super) fn rates_by_partition(partitions: &[ReconstructedPartition]) -> Vec<(usize, f64)> {
      partitions
        .iter()
        .map(|partition| {
          let gtr = partition.partition.gtr().expect("a marginal partition has a model");
          (partition.layout.columns.len(), gtr.mu())
        })
        .collect()
    }
  }
}
