#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{F81Params, GtrModel, GtrModelName};
  use crate::partition::resolve::{
    PartitionEntry, PartitionPlan, ResolvedLayouts, ResolvedPartitions, resolve_partition_layouts, resolve_partitions,
  };
  use eyre::Report;
  use helpers::{entry, layout, named};
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;

  const NO_ENTRIES: &[PartitionEntry] = &[];

  #[test]
  fn test_resolve_partitions_without_entries_gives_one_partition_per_alignment() -> Result<(), Report> {
    let model = GtrModel::from(GtrModelName::JC69);
    assert_eq!(
      ResolvedPartitions {
        plans: vec![
          PartitionPlan {
            layout: layout(1, None, 1, &[0, 1, 2]),
            model: model.clone(),
          },
          PartitionPlan {
            layout: layout(2, None, 2, &[0, 1]),
            model: model.clone(),
          },
        ],
        duplicate_names: vec![],
      },
      resolve_partitions(NO_ENTRIES, &[3, 2], &model, Some(AlphabetName::Nuc))?
    );
    Ok(())
  }

  #[test]
  fn test_resolve_partitions_codon_positions_with_a_remainder() -> Result<(), Report> {
    let hky85 = GtrModel::from(GtrModelName::HKY85);
    let entries = [
      named("cp1", entry(None, Some(&["4-12\\3"]), GtrModel::Infer)),
      named("cp2", entry(None, Some(&["5-12\\3"]), GtrModel::Infer)),
      named("cp3", entry(None, Some(&["6-12\\3"]), hky85.clone())),
      named("noncoding", entry(None, None, GtrModel::Infer)),
    ];
    assert_eq!(
      ResolvedPartitions {
        plans: vec![
          PartitionPlan {
            layout: layout(1, Some("cp1"), 1, &[3, 6, 9]),
            model: GtrModel::Infer,
          },
          PartitionPlan {
            layout: layout(2, Some("cp2"), 1, &[4, 7, 10]),
            model: GtrModel::Infer,
          },
          PartitionPlan {
            layout: layout(3, Some("cp3"), 1, &[5, 8, 11]),
            model: hky85,
          },
          PartitionPlan {
            layout: layout(4, Some("noncoding"), 1, &[0, 1, 2]),
            model: GtrModel::Infer,
          },
        ],
        duplicate_names: vec![],
      },
      resolve_partitions(
        &entries,
        &[12],
        &GtrModel::from(GtrModelName::JC69),
        Some(AlphabetName::Nuc)
      )?
    );
    Ok(())
  }

  #[test]
  fn test_resolve_partition_layouts_sort_columns_and_take_the_remainder_per_alignment() -> Result<(), Report> {
    let entries = [
      entry(Some(2), Some(&["7-9", "1-2"]), GtrModel::Infer),
      entry(Some(1), None, GtrModel::Infer),
      entry(Some(2), None, GtrModel::Infer),
    ];
    assert_eq!(
      ResolvedLayouts {
        layouts: vec![
          layout(1, None, 2, &[0, 1, 6, 7, 8]),
          layout(2, None, 1, &[0, 1, 2, 3]),
          layout(3, None, 2, &[2, 3, 4, 5]),
        ],
        duplicate_names: vec![],
      },
      resolve_partition_layouts(&entries, &[4, 9])?
    );
    Ok(())
  }

  #[test]
  fn test_resolve_partition_layouts_return_the_repeated_names() -> Result<(), Report> {
    let entries = [
      named("genes", entry(None, Some(&["1-3"]), GtrModel::Infer)),
      named("rest", entry(None, Some(&["4"]), GtrModel::Infer)),
      named("genes", entry(None, Some(&["5-6"]), GtrModel::Infer)),
      named("genes", entry(None, None, GtrModel::Infer)),
    ];
    assert_eq!(
      vec!["genes"],
      resolve_partition_layouts(&entries, &[8])?.duplicate_names
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_alignment_file(   &[],     vec![entry(None,    Some(&["1-3"]),   GtrModel::Infer)],                                                   "--partition needs --alignment: each entry takes sites of an --alignment file, but the run has no --alignment file")]
  #[case::alignment_too_large( &[9, 9], vec![named("a", entry(Some(3), None, GtrModel::Infer)), entry(Some(1), None, GtrModel::Infer)],        "partition 1 'a' names alignment 3, but --alignment lists 2 files")]
  #[case::alignment_zero(      &[9],    vec![entry(Some(0), None,             GtrModel::Infer)],                                                   "partition 1 names alignment 0, but --alignment lists 1 file")]
  #[case::alignment_required(  &[9, 9], vec![entry(Some(1), None, GtrModel::Infer), named("b", entry(None, None, GtrModel::Infer))],        "partition 2 'b' needs `alignment`, the number of its --alignment file, because --alignment lists 2 files")]
  #[case::alignment_uncovered( &[9, 9], vec![entry(Some(1), None,             GtrModel::Infer)],                                                   "alignment 2 has no partition entry: with --partition, each --alignment file needs at least one entry")]
  #[case::range_outside(       &[9],    vec![entry(None,    Some(&["5-10"]),  GtrModel::Infer), entry(None, None, GtrModel::Infer)],             "partition 1: the site range '5-10' ends after the last site 9 of alignment 1")]
  #[case::site_in_two_entries( &[9],    vec![named("cp1", entry(None, Some(&["1-9\\3"]), GtrModel::Infer)), named("x", entry(None, Some(&["4-6"]), GtrModel::Infer))], "partition 2 'x' takes site 4, which partition 1 'cp1' already takes")]
  #[case::site_twice_in_entry( &[9],    vec![entry(None,    Some(&["1-5", "5-9"]), GtrModel::Infer)],                                              "partition 1 lists site 5 twice")]
  #[case::two_remainders(      &[9],    vec![entry(None,    Some(&["1-3"]),   GtrModel::Infer), named("x", entry(None, None, GtrModel::Infer)), named("y", entry(None, None, GtrModel::Infer))], "partition 2 'x' and partition 3 'y' both leave out `sites` in alignment 1, but only one entry of an alignment can take its remaining sites")]
  #[case::uncovered_sites(     &[9],    vec![entry(None,    Some(&["1-9\\3", "2-9\\3"]), GtrModel::Infer)],                                       "no partition takes 3 sites of alignment 1, the first is site 3: list them in an entry, or add an entry without `sites`, which takes the remaining sites")]
  #[case::one_uncovered_site(  &[9],    vec![entry(None,    Some(&["1-8"]),   GtrModel::Infer)],                                                   "no partition takes 1 site of alignment 1, the first is site 9: list them in an entry, or add an entry without `sites`, which takes the remaining sites")]
  #[case::empty_sites(         &[9],    vec![entry(None,    Some(&[]),        GtrModel::Infer), entry(None, None, GtrModel::Infer)],             "partition 1 has an empty `sites` list: list its sites, or leave out `sites` to take the remaining sites of its alignment")]
  #[case::empty_remainder(     &[9],    vec![entry(None,    Some(&["1-9"]),   GtrModel::Infer), named("rest", entry(None, None, GtrModel::Infer))], "partition 2 'rest' leaves out `sites`, but the other entries take every site of alignment 1")]
  #[trace]
  fn test_resolve_partition_layouts_errors(
    #[case] lengths: &[usize],
    #[case] entries: Vec<PartitionEntry>,
    #[case] expected: &str,
  ) {
    assert_error!(resolve_partition_layouts(&entries, lengths), expected);
  }

  #[test]
  fn test_resolve_partition_layouts_reject_an_alignment_without_sites() {
    assert_error!(
      resolve_partition_layouts(NO_ENTRIES, &[4, 0]),
      "alignment 2 has no sites"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k80_aa(    GtrModelName::K80,   AlphabetName::Aa,       "partition 1 has the model k80, which needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa")]
  #[case::hky85_aa(  GtrModelName::HKY85, AlphabetName::Aa,       "partition 1 has the model hky85, which needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa")]
  #[case::t92_aa(    GtrModelName::T92,   AlphabetName::AaNoStop, "partition 1 has the model t92, which needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa-no-stop")]
  #[case::tn93_aa(   GtrModelName::TN93,  AlphabetName::Aa,       "partition 1 has the model tn93, which needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa")]
  #[case::jtt92_aa(  GtrModelName::Jtt92, AlphabetName::Aa,       "partition 1 has the model jtt92, which needs the alphabet aa-no-stop because its matrix holds the 20 amino acids without the stop '*', but the alphabet is aa")]
  #[case::jtt92_nuc( GtrModelName::Jtt92, AlphabetName::Nuc,      "partition 1 has the model jtt92, which needs the alphabet aa-no-stop because its matrix holds the 20 amino acids without the stop '*', but the alphabet is nuc")]
  #[trace]
  fn test_resolve_partitions_reject_a_model_without_its_alphabet_in_runs_without_entries(
    #[case] model: GtrModelName,
    #[case] alphabet: AlphabetName,
    #[case] expected: &str,
  ) {
    assert_error!(resolve_partitions(NO_ENTRIES, &[30], &GtrModel::from(model), Some(alphabet)), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::tn93_aa(       GtrModelName::TN93,  Some(AlphabetName::Aa), "the model tn93 needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa")]
  #[case::jtt92_no_stop( GtrModelName::Jtt92, Some(AlphabetName::Nuc), "the model jtt92 needs the alphabet aa-no-stop because its matrix holds the 20 amino acids without the stop '*', but the alphabet is nuc")]
  #[trace]
  fn test_resolve_partitions_check_the_model_of_a_run_without_alignments(
    #[case] model: GtrModelName,
    #[case] alphabet: Option<AlphabetName>,
    #[case] expected: &str,
  ) {
    assert_error!(resolve_partitions(NO_ENTRIES, &[], &GtrModel::from(model), alphabet), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::jtt92_own_alphabet(GtrModelName::Jtt92, None                  )]
  #[case::hky85_own_alphabet(GtrModelName::HKY85, None                  )]
  #[case::jc69_aa(           GtrModelName::JC69,  Some(AlphabetName::Aa))]
  #[case::infer_aa(          GtrModelName::Infer, Some(AlphabetName::Aa))]
  #[trace]
  fn test_resolve_partitions_accept_the_model_of_a_run_without_alignments(
    #[case] model: GtrModelName,
    #[case] alphabet: Option<AlphabetName>,
  ) -> Result<(), Report> {
    assert_eq!(
      ResolvedPartitions { plans: vec![], duplicate_names: vec![] },
      resolve_partitions(NO_ENTRIES, &[], &GtrModel::from(model), alphabet)?
    );
    Ok(())
  }

  #[test]
  fn test_resolve_partitions_name_the_entry_of_a_model_without_its_alphabet() {
    let entries = [
      named("cp12", entry(None, Some(&["1-30\\3", "2-30\\3"]), GtrModel::Infer)),
      named("cp3", entry(None, None, GtrModel::from(GtrModelName::HKY85))),
    ];
    assert_error!(
      resolve_partitions(&entries, &[30], &GtrModel::Infer, Some(AlphabetName::Aa)),
      "partition 2 'cp3' has the model hky85, which needs the alphabet nuc because it sets rates between nucleotides, \
       but the alphabet is aa"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::infer_aa(     GtrModelName::Infer, AlphabetName::Aa,       None    )]
  #[case::jc69_aa(      GtrModelName::JC69,  AlphabetName::Aa,       Some(21))]
  #[case::f81_aa(       GtrModelName::F81,   AlphabetName::Aa,       Some(21))]
  #[case::jc69_nuc(     GtrModelName::JC69,  AlphabetName::Nuc,      Some(4) )]
  #[case::hky85_nuc(    GtrModelName::HKY85, AlphabetName::Nuc,      Some(4) )]
  #[case::jtt92_no_stop(GtrModelName::Jtt92, AlphabetName::AaNoStop, Some(20))]
  #[trace]
  fn test_resolve_partitions_accept_a_model_that_fits_the_alphabet(
    #[case] model: GtrModelName,
    #[case] alphabet: AlphabetName,
    #[case] n_states: Option<usize>,
  ) -> Result<(), Report> {
    let ResolvedPartitions { plans, .. } = resolve_partitions(NO_ENTRIES, &[30], &GtrModel::from(model), Some(alphabet))?;
    let built = plans
      .iter()
      .filter(|plan| !plan.model.is_infer())
      .map(|plan| Ok(plan.model.build(alphabet)?.pi().len()))
      .collect::<Result<Vec<_>, Report>>()?;
    assert_eq!(n_states.into_iter().collect::<Vec<_>>(), built);
    Ok(())
  }

  #[test]
  fn test_resolve_partitions_reject_frequencies_of_another_alphabet() {
    let f81 = GtrModel::F81(F81Params {
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      ..F81Params::default()
    });
    let entries = [named("protein", entry(None, None, f81))];
    assert_error!(
      resolve_partitions(&entries, &[30], &GtrModel::Infer, Some(AlphabetName::Aa)),
      "partition 1 'protein' has the model f81 with 4 values in pi, but the alphabet aa has 21 states, so pi needs 21"
    );
  }

  mod helpers {
    use crate::gtr::get_gtr::GtrModel;
    use crate::partition::resolve::{PartitionEntry, PartitionLayout};
    use crate::partition::sites::SiteRange;

    pub(super) fn entry(alignment: Option<usize>, sites: Option<&[&str]>, model: GtrModel) -> PartitionEntry {
      PartitionEntry {
        name: None,
        alignment,
        sites: sites.map(|sites| sites.iter().map(|text| text.parse::<SiteRange>().unwrap()).collect()),
        model,
      }
    }

    pub(super) fn named(name: &str, entry: PartitionEntry) -> PartitionEntry {
      PartitionEntry {
        name: Some(name.to_owned()),
        ..entry
      }
    }

    pub(super) fn layout(number: usize, name: Option<&str>, alignment: usize, columns: &[usize]) -> PartitionLayout {
      PartitionLayout {
        number,
        name: name.map(str::to_owned),
        alignment,
        columns: columns.to_vec(),
      }
    }
  }
}
