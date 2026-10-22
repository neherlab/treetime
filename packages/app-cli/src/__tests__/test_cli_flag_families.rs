#[cfg(test)]
pub(crate) mod tests {
  use crate::__tests__::test_cli_support::tests::helpers::{
    AN, ANALYSIS, CL, COMPLETIONS, EVERY, EXAMPLES_GET, HELP_MARKDOWN, HO, MU, OP, PI, PR, RUN, SCHEMA_WRITE, TT,
    find_command, leaf_commands, long_flags, treetime_command,
  };
  use helpers::{expected_flags, per_file_output_flags};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;

  pub(crate) struct Family {
    pub(crate) name: &'static str,
    pub(crate) flags: &'static [&'static str],
    pub(crate) commands: &'static [&'static str],
  }

  #[rustfmt::skip]
  pub(crate) const FAMILIES: &[Family] = &[
    Family { name: "Global",                   commands: EVERY,                flags: &["--verbose", "--quiet", "--verbosity", "--silent", "--no-progress", "--color", "--help"] },
    Family { name: "Run",                      commands: RUN,                  flags: &["--config", "--jobs", "--output-all"] },
    Family { name: "Preview",                  commands: &[TT, AN, CL, MU, HO, OP, PR, PI, EXAMPLES_GET, SCHEMA_WRITE],
                                                                               flags: &["--dry-run"] },
    Family { name: "Output",                   commands: ANALYSIS,             flags: &["--output-selection"] },
    Family { name: "Tree input and output",    commands: ANALYSIS,             flags: &["--tree", "--newick-dialect", "--newick-annotations", "--topology-order",
                                                                                        "--topology-order-tree", "--topology-order-list", "--topology-order-aggregate"] },
    Family { name: "Sequence input",           commands: &[TT, AN, HO, OP, PR], flags: &["--alignment", "--alphabet", "--gap-fill", "--allow-missing-sequences", "--zero-based"] },
    Family { name: "Substitution model",       commands: &[TT, AN, HO, OP],    flags: &["--model", "--model-params", "--representation"] },
    Family { name: "Reconstruction method",    commands: &[AN, HO],            flags: &["--method-anc", "--sample-from-profile"] },
    Family { name: "Tip sequences",            commands: &[TT, AN],            flags: &["--include-leaves", "--reconstruct-tip-states", "--report-ambiguous"] },
    Family { name: "Tip imputation",           commands: &[TT, AN, HO],        flags: &["--impute-missing-data"] },
    Family { name: "Ambiguity report",         commands: &[TT, AN, MU],        flags: &["--ambiguity-threshold", "--ambiguity-outputs"] },
    Family { name: "Branch-length estimation", commands: &[TT, OP],            flags: &["--max-iter", "--no-indels", "--partition-rates", "--divergence-units"] },
    Family { name: "Metadata",                 commands: &[TT, CL, MU],        flags: &["--metadata", "--metadata-id-columns", "--metadata-delimiters"] },
    Family { name: "Dates",                    commands: &[TT, CL],            flags: &["--date-column", "--date-format"] },
    Family { name: "Rooting",                  commands: &[TT, CL, OP],        flags: &["--reroot", "--outgroup", "--keep-root"] },
    Family { name: "Clock filter",             commands: &[TT, CL],            flags: &["--clock-filter", "--allow-negative-rate"] },
    Family { name: "Random steps",             commands: &[TT, AN, HO],        flags: &["--seed"] },
    Family { name: "Uncertainty",              commands: &[TT, MU],            flags: &["--confidence"] },
  ];

  #[rustfmt::skip]
  const OWN_FLAGS: &[(&str, &[&str])] = &[
    (TT,            &["--sequence-length", "--clock-rate", "--clock-std-dev", "--relax", "--branch-length-mode",
                      "--keep-polytomies", "--max-grid-points", "--coalescent", "--n-skyline", "--skyline-stiffness",
                      "--coalescent-confidence", "--gen-per-year"]),
    (AN,            &["--translations", "--cds-selection", "--annotation", "--aa-root-sequence", "--aa-model"]),
    (CL,            &[]),
    (MU,            &["--attribute", "--weights", "--missing-data", "--pc", "--missing-weights-threshold",
                      "--model-iterations", "--sampling-bias-correction", "--smooth-initial-pi",
                      "--filter-uninformative-root"]),
    (HO,            &["--const", "--rescale", "--detailed", "--drms", "--num-mut"]),
    (OP,            &["--tolerance", "--damping", "--initial-branch-lengths", "--optimizer", "--no-collapse-zero-branches",
                      "--no-merge-shared-mutations", "--no-hoist-reversions"]),
    (PR,            &["--collapse-short", "--collapse-empty", "--merge-shared-mutations", "--remove", "--remove-file"]),
    (PI,            &["--steps"]),
    (EXAMPLES_GET,  &["--output-dir", "--url"]),
    (SCHEMA_WRITE,  &["--for", "--output"]),
    (COMPLETIONS,   &[]),
    (HELP_MARKDOWN, &[]),
  ];

  #[test]
  fn test_cli_flag_families_cover_every_command() {
    let root = treetime_command();
    let expected: BTreeSet<String> = EVERY.iter().map(|path| (*path).to_owned()).collect();

    let actual: BTreeSet<String> = leaf_commands(&root).into_iter().map(|(path, _)| path).collect();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_cli_flag_families_list_each_flag_once_per_command() {
    let expected: Vec<(&str, &str)> = vec![];

    let actual = EVERY
      .iter()
      .flat_map(|command| {
        FAMILIES
          .iter()
          .filter(|family| family.commands.contains(command))
          .flat_map(|family| family.flags.iter())
          .chain(
            OWN_FLAGS
              .iter()
              .filter(move |(path, _)| path == command)
              .flat_map(|(_, flags)| flags.iter()),
          )
          .duplicates()
          .map(|flag| (*command, *flag))
      })
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree(     TT)]
  #[case::ancestral(    AN)]
  #[case::clock(        CL)]
  #[case::mugration(    MU)]
  #[case::homoplasy(    HO)]
  #[case::optimize(     OP)]
  #[case::prune(        PR)]
  #[case::pipeline(     PI)]
  #[case::examples_get( EXAMPLES_GET)]
  #[case::schema_write( SCHEMA_WRITE)]
  #[case::completions(  COMPLETIONS)]
  #[case::help_markdown(HELP_MARKDOWN)]
  #[trace]
  fn test_cli_flag_families_give_each_command_exactly_its_flags(#[case] path: &str) {
    let root = treetime_command();
    let command = find_command(&root, path).unwrap();
    let expected = expected_flags(path, FAMILIES, OWN_FLAGS);

    let per_file_outputs = per_file_output_flags(command);
    let actual: BTreeSet<String> = long_flags(command)
      .into_iter()
      .filter(|flag| !(ANALYSIS.contains(&path) && per_file_outputs.contains(flag)))
      .collect();

    assert_eq!(expected, actual);
  }

  pub(crate) mod helpers {
    use super::Family;
    use clap::Command;
    use std::any::TypeId;
    use std::collections::BTreeSet;
    use std::path::PathBuf;

    pub(super) fn expected_flags(path: &str, families: &[Family], own: &[(&str, &[&str])]) -> BTreeSet<String> {
      families
        .iter()
        .filter(|family| family.commands.contains(&path))
        .flat_map(|family| family.flags.iter())
        .chain(
          own
            .iter()
            .filter(|(command, _)| *command == path)
            .flat_map(|(_, flags)| flags.iter()),
        )
        .map(|flag| (*flag).to_owned())
        .collect()
    }

    pub(crate) fn per_file_output_flags(command: &Command) -> BTreeSet<String> {
      command
        .get_arguments()
        .filter(|arg| arg.get_value_parser().type_id() == TypeId::of::<PathBuf>())
        .filter_map(|arg| arg.get_long())
        .filter(|long| long.starts_with("output-") && *long != "output-all")
        .map(|long| format!("--{long}"))
        .collect()
    }
  }
}
