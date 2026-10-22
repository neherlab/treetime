#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::RecordingLog;
  use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
  use crate::commands::shared::argument_error::ArgumentError;
  use crate::commands::shared::partition::{partition_entries, resolve_run_partitions};
  use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
  use crate::yaml::yaml_read_str;
  use clap::Parser;
  use clap::error::ErrorKind;
  use eyre::Report;
  use helpers::{argument_error, sites};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::json;
  use treetime::alphabet::alphabet::AlphabetName;
  use treetime::gtr::get_gtr::{GtrModel, GtrModelName, HKY85Params};
  use treetime::partition::resolve::{PartitionEntry, PartitionPlan};
  use treetime::progress::{RunWarning, RunWarningKind};
  use treetime_utils::io::json::to_json_value;
  use treetime_utils::{o, vec_of_owned};

  const ANCESTRAL: [&str; 5] = ["ancestral", "--tree", "t.nwk", "--alignment", "a.fasta"];

  #[test]
  fn test_partition_entries_carry_their_models() -> Result<(), Report> {
    let raw = TreetimeAncestralArgsRaw::try_parse_from(
      [
        &ANCESTRAL[..],
        &[
          "--partition",
          "{name: cp12, sites: ['1-1500\\3', '2-1500\\3']}",
          "--partition",
          "{name: cp3, sites: ['3-1500\\3'], model: hky85, model_params: ['kappa=4']}",
          "--partition",
          r#"{"alignment": 1, "model": "jc69"}"#,
        ],
      ]
      .concat(),
    )?;
    assert_eq!(
      vec![
        PartitionEntry {
          name: Some(o!("cp12")),
          alignment: None,
          sites: Some(sites(&["1-1500\\3", "2-1500\\3"])),
          model: GtrModel::Infer,
        },
        PartitionEntry {
          name: Some(o!("cp3")),
          alignment: None,
          sites: Some(sites(&["3-1500\\3"])),
          model: GtrModel::HKY85(HKY85Params {
            kappa: 4.0,
            ..HKY85Params::default()
          }),
        },
        PartitionEntry {
          name: None,
          alignment: Some(1),
          sites: None,
          model: GtrModel::from(GtrModelName::JC69),
        },
      ],
      TreetimeAncestralArgs::try_from(raw)?.sequences.partition_entries
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::params_on_infer(       "{name: cp3, model_params: ['kappa=4']}",             "partition 1 'cp3': model_params needs a named model: infer fits the model from the data and takes no parameters, but got 'kappa=4'")]
  #[case::params_without_value(  "{model: hky85, model_params: ['kappa']}",            "partition 1: model_params takes key=value items, but got 'kappa'")]
  #[case::param_of_another_model("{name: x, model: jc69, model_params: ['kappa=2']}",  "partition 1 'x': model_params item 'kappa=2': the model jc69 has no parameter 'kappa'; its parameters are mu")]
  #[case::param_twice(           "{model: k80, model_params: ['kappa=2', 'kappa=3']}", "partition 1: model_params sets kappa twice: 'kappa=2' and 'kappa=3'")]
  #[case::negative_rate(         "{model: k80, model_params: ['mu=-1']}",              "partition 1: model_params item 'mu=-1': mu must be a non-negative finite number")]
  #[trace]
  fn test_partition_entry_errors_name_the_entry(#[case] entry: &str, #[case] expected: &str) {
    let raw = TreetimeAncestralArgsRaw::try_parse_from([&ANCESTRAL[..], &["--partition", entry]].concat()).unwrap();
    assert_eq!(
      Some((o!(expected), vec_of_owned!["partition"])),
      argument_error(TreetimeAncestralArgs::try_from(raw))
    );
  }

  #[test]
  fn test_partition_entry_errors_name_their_number() {
    let raw = TreetimeAncestralArgsRaw::try_parse_from(
      [
        &ANCESTRAL[..],
        &[
          "--partition",
          "{sites: ['1-9']}",
          "--partition",
          "{name: rest, model_params: ['mu=2']}",
        ],
      ]
      .concat(),
    )
    .unwrap();
    assert_eq!(
      Some((
        o!(
          "partition 2 'rest': model_params needs a named model: infer fits the model from the data and takes no \
           parameters, but got 'mu=2'"
        ),
        vec_of_owned!["partition"]
      )),
      argument_error(TreetimeAncestralArgs::try_from(raw))
    );
  }

  #[test]
  fn test_partition_flags_and_config_give_the_same_entries() -> Result<(), Report> {
    let flags = TreetimeTimetreeArgsRaw::try_parse_from([
      "timetree",
      "--partition",
      "{name: cp1, sites: ['101-1600\\3']}",
      "--partition",
      "{name: cp3, sites: ['103-1600\\3', 1-100], model: hky85, model_params: ['kappa=4']}",
      "--partition",
      "{name: noncoding}",
    ])?;
    let config: TreetimeTimetreeArgsRaw = yaml_read_str(indoc! {r"
      partition:
        - {name: cp1, sites: ['101-1600\3']}
        - name: cp3
          sites:
            - '103-1600\3'
            - 1-100
          model: hky85
          model_params: ['kappa=4']
        - {name: noncoding}
    "})?;
    assert_eq!(flags.partition.partition, config.partition.partition);
    assert_eq!(3, flags.partition.partition.len());
    Ok(())
  }

  #[test]
  fn test_partition_entries_serialize_as_written() -> Result<(), Report> {
    let config: TreetimeTimetreeArgsRaw = yaml_read_str(indoc! {r"
      partition:
        - {name: cp3, sites: ['103-1600\3', '5'], model: hky85, model_params: ['kappa=4']}
        - {alignment: 2}
        - {}
    "})?;
    assert_eq!(
      json!([
        {"name": "cp3", "sites": ["103-1600\\3", "5"], "model": "hky85", "model_params": ["kappa=4"]},
        {"alignment": 2},
        {},
      ]),
      to_json_value(&config.partition.partition)?
    );
    Ok(())
  }

  #[test]
  fn test_partition_flag_rejects_an_invalid_site_range() {
    let error = TreetimeTimetreeArgsRaw::try_parse_from(["timetree", "--partition", "{sites: ['9-5']}"]).unwrap_err();
    assert_eq!(
      (
        ErrorKind::ValueValidation,
        o!(
          "error: invalid value '{sites: ['9-5']}' for '--partition <ENTRY>': Unexpected: invalid value: site range \
           '9-5': the first site 9 is after the last site 5 at line 1 column 10 (path: sites[0])\n\nFor more \
           information, try '--help'.\n"
        )
      ),
      (error.kind(), error.to_string())
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::model(       &["--model", "jc69"])]
  #[case::model_params(&["--model-params", "mu=2"])]
  #[trace]
  fn test_partition_conflicts_with_the_model_flags(#[case] model: &[&str]) {
    let error = TreetimeAncestralArgsRaw::try_parse_from([&ANCESTRAL[..], &["--partition", "{}"], model].concat())
      .unwrap_err();
    assert_eq!(ErrorKind::ArgumentConflict, error.kind());
  }

  #[test]
  fn test_partition_repeated_names_raise_one_warning() -> Result<(), Report> {
    let raw = TreetimeAncestralArgsRaw::try_parse_from(
      [
        &ANCESTRAL[..],
        &[
          "--partition",
          "{name: gene, sites: ['1-3']}",
          "--partition",
          "{name: gene, sites: ['4-6']}",
          "--partition",
          "{name: rest}",
        ],
      ]
      .concat(),
    )?;
    let args = TreetimeAncestralArgs::try_from(raw)?;
    let log = RecordingLog::default();
    let plans = resolve_run_partitions(
      &args.sequences.partition_entries,
      &GtrModel::Infer,
      &[9],
      Some(AlphabetName::Nuc),
      &log,
    )?;

    assert_eq!(
      (
        vec![vec![0, 1, 2], vec![3, 4, 5], vec![6, 7, 8]],
        vec![RunWarning {
          kind: RunWarningKind::DuplicatePartitionNames,
          message: o!(
            "More than one partition entry is named gene. TreeTime tells these partitions apart by their numbers."
          ),
          names: vec_of_owned!["gene"],
        }]
      ),
      (
        plans
          .into_iter()
          .map(|PartitionPlan { layout, .. }| layout.columns)
          .collect::<Vec<_>>(),
        log.into_warnings()
      )
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::layout_with_entries(&["{sites: ['1-4']}"], GtrModelName::Infer, AlphabetName::Nuc, "no partition takes 5 sites of alignment 1, the first is site 5: list them in an entry, or add an entry without `sites`, which takes the remaining sites", &["partition", "alignment"])]
  #[case::model_with_entries( &["{model: k80}"],     GtrModelName::Infer, AlphabetName::Aa,  "partition 1 has the model k80, which needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa", &["partition", "alphabet"])]
  #[case::model_of_the_flags( &[],                   GtrModelName::K80,   AlphabetName::Aa,  "partition 1 has the model k80, which needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa", &["model", "model_params", "alphabet"])]
  #[trace]
  fn test_partition_resolution_errors_are_argument_errors(
    #[case] specs: &[&str],
    #[case] model: GtrModelName,
    #[case] alphabet: AlphabetName,
    #[case] message: &str,
    #[case] settings: &[&str],
  ) -> Result<(), Report> {
    let specs = specs.iter().map(|spec| yaml_read_str(spec)).collect::<Result<Vec<_>, Report>>()?;
    let entries = partition_entries(&specs)?;
    assert_eq!(
      Some((o!(message), settings.iter().map(|setting| (*setting).to_owned()).collect())),
      argument_error(resolve_run_partitions(
        &entries,
        &GtrModel::from(model),
        &[9],
        Some(alphabet),
        &RecordingLog::default()
      ))
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::hky85_aa(      GtrModelName::HKY85, Some(AlphabetName::Aa), Some(o!("the model hky85 needs the alphabet nuc because it sets rates between nucleotides, but the alphabet is aa")))]
  #[case::jtt92_aa(      GtrModelName::Jtt92, Some(AlphabetName::Aa), Some(o!("the model jtt92 needs the alphabet aa-no-stop because its matrix holds the 20 amino acids without the stop '*', but the alphabet is aa")))]
  #[case::jtt92_default( GtrModelName::Jtt92, None,                   None)]
  #[case::hky85_default( GtrModelName::HKY85, None,                   None)]
  #[case::infer_aa(      GtrModelName::Infer, Some(AlphabetName::Aa), None)]
  #[trace]
  fn test_partition_resolution_checks_the_model_of_a_run_without_alignments(
    #[case] model: GtrModelName,
    #[case] alphabet: Option<AlphabetName>,
    #[case] message: Option<String>,
  ) {
    assert_eq!(
      message.map(|message| (message, vec_of_owned!["model", "model_params", "alphabet"])),
      argument_error(resolve_run_partitions(
        &[],
        &GtrModel::from(model),
        &[],
        alphabet,
        &RecordingLog::default()
      ))
    );
  }

  mod helpers {
    use super::ArgumentError;
    use eyre::Report;
    use treetime::partition::sites::SiteRange;

    pub(super) fn argument_error<T>(result: Result<T, Report>) -> Option<(String, Vec<String>)> {
      result.err().map(|report| {
        let error = report
          .downcast_ref::<ArgumentError>()
          .unwrap_or_else(|| panic!("expected an argument error, got: {report:?}"));
        (error.message.clone(), error.settings.clone())
      })
    }

    pub(super) fn sites(texts: &[&str]) -> Vec<SiteRange> {
      texts.iter().map(|text| text.parse().unwrap()).collect()
    }
  }
}
