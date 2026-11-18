#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{RecordingLog, argument_error};
  use crate::commands::shared::partition_by::{
    AnnotatedAlignment, RateCategoryRule, codon_position_entries, gene_entries, k_means_thresholds,
    rate_category_buckets, rate_category_entries,
  };
  use eyre::Report;
  use helpers::{
    annotation, bucket_sites, cds, entry_sites, exhaustive_min_squared_deviations, region, squared_deviations_of_split,
  };
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime::gtr::get_gtr::GtrModel;
  use treetime::progress::{RunWarning, RunWarningKind};
  use treetime_io::gff::{GffPhase, GffStrand};
  use treetime_utils::{o, vec_of_owned};

  #[test]
  fn test_codon_positions_of_forward_and_reverse_cds() -> Result<(), Report> {
    let annotation = annotation(
      vec![
        cds("A", GffStrand::Forward, &[(1, 9, GffPhase::Zero)]),
        cds("B", GffStrand::Reverse, &[(11, 19, GffPhase::Zero)]),
      ],
      20,
    );
    let actual = codon_position_entries(
      &[AnnotatedAlignment {
        length: 20,
        annotation: &annotation,
      }],
      &GtrModel::Infer,
      &RecordingLog::default(),
    )?;

    let expected = vec![
      (o!("cp1"), Some(vec![1, 4, 7, 13, 16, 19])),
      (o!("cp2"), Some(vec![2, 5, 8, 12, 15, 18])),
      (o!("cp3"), Some(vec![3, 6, 9, 11, 14, 17])),
      (o!("noncoding"), None),
    ];
    assert_eq!(
      expected,
      entry_sites(&actual),
      "codon positions counted from the 5' end of each CDS"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::phase_one(&[(1, 5, GffPhase::Zero), (10, 16, GffPhase::One)], vec![1, 4, 11, 14], vec![2, 5, 12, 15], vec![3, 10, 13, 16])]
  #[case::phase_two(&[(1, 4, GffPhase::Zero), (10, 17, GffPhase::Two)], vec![1, 4, 12, 15], vec![2, 10, 13, 16], vec![3, 11, 14, 17])]
  #[trace]
  fn test_codon_positions_follow_the_phase_of_each_segment(
    #[case] segments: &[(usize, usize, GffPhase)],
    #[case] cp1: Vec<usize>,
    #[case] cp2: Vec<usize>,
    #[case] cp3: Vec<usize>,
  ) -> Result<(), Report> {
    let annotation = annotation(vec![cds("A", GffStrand::Forward, segments)], 20);
    let actual = codon_position_entries(
      &[AnnotatedAlignment { length: 20, annotation: &annotation }],
      &GtrModel::Infer,
      &RecordingLog::default(),
    )?;

    let expected = vec![
      (o!("cp1"), Some(cp1)),
      (o!("cp2"), Some(cp2)),
      (o!("cp3"), Some(cp3)),
      (o!("noncoding"), None),
    ];
    assert_eq!(expected, entry_sites(&actual), "the phase of a GFF3 row counts the bases before its first codon");
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::longer_wins(   vec![cds("A", GffStrand::Forward, &[(1, 12, GffPhase::Zero)]), cds("B", GffStrand::Forward, &[(10, 15, GffPhase::Zero)])], vec![(o!("A"), Some((1..=12).collect())), (o!("B"), Some(vec![13, 14, 15])), (o!("noncoding"), None)])]
  #[case::lower_start(   vec![cds("B", GffStrand::Forward, &[(4, 9, GffPhase::Zero)]),  cds("A", GffStrand::Forward, &[(1, 6, GffPhase::Zero)])],   vec![(o!("A"), Some((1..=6).collect())),  (o!("B"), Some(vec![7, 8, 9])),    (o!("noncoding"), None)])]
  #[case::label_breaks(  vec![cds("Y", GffStrand::Forward, &[(1, 6, GffPhase::Zero)]),  cds("X", GffStrand::Forward, &[(4, 9, GffPhase::Zero)]), cds("W", GffStrand::Forward, &[(7, 12, GffPhase::Zero)]), cds("V", GffStrand::Forward, &[(7, 12, GffPhase::Zero)])], vec![(o!("Y"), Some((1..=6).collect())), (o!("X"), Some(vec![7, 8, 9])), (o!("V"), Some(vec![10, 11, 12])), (o!("noncoding"), None)])]
  #[trace]
  fn test_genes_give_shared_sites_to_the_highest_ranked_cds(
    #[case] features: Vec<treetime_io::gff::GffCdsFeature>,
    #[case] expected: Vec<(String, Option<Vec<usize>>)>,
  ) -> Result<(), Report> {
    let annotation = annotation(features, 20);
    let actual = gene_entries(
      &[AnnotatedAlignment { length: 20, annotation: &annotation }],
      &GtrModel::Infer,
      &RecordingLog::default(),
    )?;
    assert_eq!(expected, entry_sites(&actual), "longest CDS first, then the lowest start, then the label");
    Ok(())
  }

  #[test]
  fn test_genes_warn_about_a_cds_inside_a_longer_cds() -> Result<(), Report> {
    let annotation = annotation(
      vec![
        cds("N", GffStrand::Forward, &[(1, 12, GffPhase::Zero)]),
        cds("ORF9b", GffStrand::Forward, &[(4, 9, GffPhase::Zero)]),
      ],
      15,
    );
    let log = RecordingLog::default();
    let actual = gene_entries(
      &[AnnotatedAlignment {
        length: 15,
        annotation: &annotation,
      }],
      &GtrModel::Infer,
      &log,
    )?;

    assert_eq!(
      vec![(o!("N"), Some((1..=12).collect())), (o!("noncoding"), None)],
      entry_sites(&actual)
    );
    assert_eq!(
      vec![RunWarning {
        kind: RunWarningKind::CdsWithoutPartitionSites,
        message: o!(
          "--partition-by makes no partition for CDS ORF9b: each of their sites is also covered by a longer CDS, and \
           the longest CDS that covers a site owns it."
        ),
        names: vec_of_owned!["ORF9b"],
      }],
      log.into_warnings()
    );
    Ok(())
  }

  #[test]
  fn test_codon_positions_without_unowned_sites_have_no_noncoding_partition() -> Result<(), Report> {
    let annotation = annotation(vec![cds("A", GffStrand::Forward, &[(1, 12, GffPhase::Zero)])], 12);
    let actual = codon_position_entries(
      &[AnnotatedAlignment {
        length: 12,
        annotation: &annotation,
      }],
      &GtrModel::Infer,
      &RecordingLog::default(),
    )?;
    let expected = vec![
      (o!("cp1"), Some(vec![1, 4, 7, 10])),
      (o!("cp2"), Some(vec![2, 5, 8, 11])),
      (o!("cp3"), Some(vec![3, 6, 9, 12])),
    ];
    assert_eq!(expected, entry_sites(&actual));
    Ok(())
  }

  #[test]
  fn test_genes_of_two_alignments_carry_their_alignment_numbers() -> Result<(), Report> {
    let first = annotation(vec![cds("A", GffStrand::Forward, &[(1, 6, GffPhase::Zero)])], 9);
    let second = annotation(vec![cds("A", GffStrand::Forward, &[(4, 9, GffPhase::Zero)])], 9);
    let actual = gene_entries(
      &[
        AnnotatedAlignment {
          length: 9,
          annotation: &first,
        },
        AnnotatedAlignment {
          length: 9,
          annotation: &second,
        },
      ],
      &GtrModel::Infer,
      &RecordingLog::default(),
    )?;

    let summary = actual
      .iter()
      .map(|entry| (entry.name.clone(), entry.alignment))
      .collect::<Vec<_>>();
    let expected = vec![
      (Some(o!("1/A")), Some(1)),
      (Some(o!("1/noncoding")), Some(1)),
      (Some(o!("2/A")), Some(2)),
      (Some(o!("2/noncoding")), Some(2)),
    ];
    assert_eq!(expected, summary);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::two_seqids(       helpers::two_seqids(),                      "The annotation of alignment 1 places its CDSs on the sequences 'seq1', 'seq2', but one annotation must describe the positions of one alignment, so all its CDSs need the same sequence ID")]
  #[case::no_region(        helpers::without_region(),                  "The annotation of alignment 1 has no ##sequence-region line for the sequence 'seq1': TreeTime needs '##sequence-region seq1 1 20' to check that the annotation positions are the 20 columns of the alignment")]
  #[case::region_length(    helpers::with_region(region("seq1", 1, 25)), "The annotation of alignment 1 gives the sequence 'seq1' the region 1-25, but the alignment has 20 columns: annotation positions are alignment columns only for an alignment on that reference, without insertion columns, so the region must run from 1 to 20")]
  #[case::cds_past_the_end( helpers::past_the_end(),                    "The annotation of alignment 1 has CDS 'A' with a segment that ends at position 24, after the last column 20 of the alignment")]
  #[trace]
  fn test_annotation_checks_are_argument_errors(
    #[case] annotation: treetime_io::gff::GffCdsAnnotation,
    #[case] message: &str,
  ) {
    let result = gene_entries(
      &[AnnotatedAlignment { length: 20, annotation: &annotation }],
      &GtrModel::Infer,
      &RecordingLog::default(),
    );
    assert_eq!(Some((o!(message), vec_of_owned!["annotation", "alignment"])), argument_error(result));
  }

  #[test]
  fn test_scheme_with_one_partition_is_an_argument_error() {
    let annotation = annotation(vec![cds("A", GffStrand::Forward, &[(1, 12, GffPhase::Zero)])], 12);
    let result = gene_entries(
      &[AnnotatedAlignment {
        length: 12,
        annotation: &annotation,
      }],
      &GtrModel::Infer,
      &RecordingLog::default(),
    );
    assert_eq!(
      Some((
        o!(
          "--partition-by makes only 1 partition ('A') over all --alignment files, but a partition scheme needs at \
           least two partitions"
        ),
        vec_of_owned!["partition_by"]
      )),
      argument_error(result)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::change_counts(     RateCategoryRule::ChangeCounts, vec![0, 1, 2, 5, 0],      vec![("slow", vec![0, 4]), ("medium", vec![1]),       ("fast", vec![2, 3])])]
  #[case::median_even(       RateCategoryRule::MedianSplit,  vec![0, 1, 2, 3, 4],      vec![("slow", vec![0]),    ("medium", vec![1, 2]),    ("fast", vec![3, 4])])]
  #[case::median_odd(        RateCategoryRule::MedianSplit,  vec![0, 1, 3, 5],         vec![("slow", vec![0]),    ("medium", vec![1, 2]),    ("fast", vec![3])])]
  #[case::median_equal(      RateCategoryRule::MedianSplit,  vec![0, 2, 2],            vec![("slow", vec![0]),    ("medium", vec![1, 2])])]
  #[case::k_means(           RateCategoryRule::KMeans,       vec![1, 1, 2, 10, 11, 30], vec![("slow", vec![0, 1, 2]), ("medium", vec![3, 4]), ("fast", vec![5])])]
  #[case::k_means_two_values(RateCategoryRule::KMeans,       vec![0, 0, 4],            vec![("slow", vec![0, 1]), ("medium", vec![2])])]
  #[trace]
  fn test_rate_category_buckets(
    #[case] rule: RateCategoryRule,
    #[case] counts: Vec<usize>,
    #[case] expected: Vec<(&str, Vec<usize>)>,
  ) {
    assert_eq!(expected, bucket_sites(&rate_category_buckets(rule, &counts)), "groups worked out by hand");
  }

  #[test]
  fn test_rate_categories_of_alignments_without_changes_are_an_argument_error() {
    let result = rate_category_entries(
      RateCategoryRule::ChangeCounts,
      &[vec![0, 0, 0]],
      &GtrModel::Infer,
      &RecordingLog::default(),
    );
    assert_eq!(
      Some((
        o!(
          "--partition-by makes only 1 partition ('slow') over all --alignment files, but a partition scheme needs \
           at least two partitions"
        ),
        vec_of_owned!["partition_by"]
      )),
      argument_error(result)
    );
  }

  #[test]
  fn test_rate_categories_of_two_alignments_without_changes_pass() -> Result<(), Report> {
    let entries = rate_category_entries(
      RateCategoryRule::ChangeCounts,
      &[vec![0, 0], vec![0]],
      &GtrModel::Infer,
      &RecordingLog::default(),
    )?;
    let names = entries
      .iter()
      .filter_map(|entry| entry.name.clone())
      .collect::<Vec<_>>();
    assert_eq!(vec_of_owned!["1/slow", "2/slow"], names);
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_k_means_thresholds_match_exhaustive_search(
      counts in proptest::collection::vec(0_usize..8, 3..8),
    ) {
      let distinct = counts.iter().collect::<std::collections::BTreeSet<_>>().len();
      prop_assume!(distinct >= 3);
      let thresholds = k_means_thresholds(&counts);
      let actual = squared_deviations_of_split(&counts, &thresholds);
      let expected = exhaustive_min_squared_deviations(&counts);
      prop_assert!((expected - actual).abs() <= 1e-9, "split {thresholds:?}: {actual} against the optimum {expected}");
    }
  }

  mod helpers {
    use crate::commands::shared::partition_by::RateBucket;
    use itertools::Itertools;
    use treetime::partition::resolve::PartitionEntry;
    use treetime_io::gff::{GffCdsAnnotation, GffCdsFeature, GffCdsSegment, GffPhase, GffSequenceRegion, GffStrand};
    use treetime_utils::o;

    pub(super) fn cds(label: &str, strand: GffStrand, segments: &[(usize, usize, GffPhase)]) -> GffCdsFeature {
      cds_on("seq1", label, strand, segments)
    }

    pub(super) fn cds_on(
      seqid: &str,
      label: &str,
      strand: GffStrand,
      segments: &[(usize, usize, GffPhase)],
    ) -> GffCdsFeature {
      GffCdsFeature {
        label: o!(label),
        id: None,
        name: o!(label),
        seqid: o!(seqid),
        strand,
        segments: segments
          .iter()
          .map(|&(start, end, phase)| GffCdsSegment { start, end, phase })
          .collect(),
      }
    }

    pub(super) fn region(seqid: &str, start: usize, end: usize) -> GffSequenceRegion {
      GffSequenceRegion {
        seqid: o!(seqid),
        start,
        end,
      }
    }

    pub(super) fn annotation(features: Vec<GffCdsFeature>, length: usize) -> GffCdsAnnotation {
      GffCdsAnnotation {
        features,
        sequence_regions: vec![region("seq1", 1, length)],
      }
    }

    pub(super) fn two_seqids() -> GffCdsAnnotation {
      GffCdsAnnotation {
        features: vec![
          cds_on("seq1", "A", GffStrand::Forward, &[(1, 6, GffPhase::Zero)]),
          cds_on("seq2", "B", GffStrand::Forward, &[(7, 12, GffPhase::Zero)]),
        ],
        sequence_regions: vec![region("seq1", 1, 20)],
      }
    }

    pub(super) fn without_region() -> GffCdsAnnotation {
      GffCdsAnnotation {
        features: vec![cds("A", GffStrand::Forward, &[(1, 6, GffPhase::Zero)])],
        sequence_regions: vec![],
      }
    }

    pub(super) fn with_region(region: GffSequenceRegion) -> GffCdsAnnotation {
      GffCdsAnnotation {
        features: vec![cds("A", GffStrand::Forward, &[(1, 6, GffPhase::Zero)])],
        sequence_regions: vec![region],
      }
    }

    pub(super) fn past_the_end() -> GffCdsAnnotation {
      GffCdsAnnotation {
        features: vec![cds("A", GffStrand::Forward, &[(1, 24, GffPhase::Zero)])],
        sequence_regions: vec![region("seq1", 1, 20)],
      }
    }

    pub(super) fn entry_sites(entries: &[PartitionEntry]) -> Vec<(String, Option<Vec<usize>>)> {
      entries
        .iter()
        .map(|entry| {
          let sites = entry.sites.as_ref().map(|ranges| {
            ranges
              .iter()
              .flat_map(|range| range.columns())
              .map(|column| column + 1)
              .sorted()
              .collect()
          });
          (entry.name.clone().unwrap_or_default(), sites)
        })
        .collect()
    }

    pub(super) fn bucket_sites(buckets: &[RateBucket]) -> Vec<(&str, Vec<usize>)> {
      buckets
        .iter()
        .map(|bucket| (bucket.name.as_str(), bucket.columns.iter().copied().collect()))
        .collect()
    }

    pub(super) fn squared_deviations_of_split(counts: &[usize], thresholds: &[usize]) -> f64 {
      let groups = counts
        .iter()
        .map(|&count| {
          (
            thresholds.iter().take_while(|&&threshold| count > threshold).count(),
            count,
          )
        })
        .into_group_map();
      groups.values().map(|group| squared_deviations(group)).sum()
    }

    pub(super) fn exhaustive_min_squared_deviations(counts: &[usize]) -> f64 {
      (0..3_usize.pow(u32::try_from(counts.len()).unwrap()))
        .map(|code| {
          let groups = counts
            .iter()
            .scan(code, |code, &count| {
              let group = *code % 3;
              *code /= 3;
              Some((group, count))
            })
            .into_group_map();
          groups.values().map(|group| squared_deviations(group)).sum::<f64>()
        })
        .fold(f64::INFINITY, f64::min)
    }

    fn squared_deviations(values: &[usize]) -> f64 {
      let values = values
        .iter()
        .map(|&value| f64::from(u32::try_from(value).unwrap()))
        .collect_vec();
      let mean = values.iter().sum::<f64>() / f64::from(u32::try_from(values.len()).unwrap());
      values.iter().map(|value| (value - mean).powi(2)).sum()
    }
  }
}
