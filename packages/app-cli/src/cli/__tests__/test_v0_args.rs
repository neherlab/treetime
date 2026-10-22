#[cfg(test)]
mod tests {
  use crate::cli::v0_args::{EVERY_COMMAND, ROOT, Then, V0_SPELLINGS, V0Rewrite, Values, rewrite_v0_args};
  use clap::Command;
  use helpers::{accepts, attached, find_command, long_flag, non_utf8, rewrite, rewritten_args, treetime_command};
  use itertools::Itertools;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::ffi::OsString;
  use std::iter;
  use treetime_utils::{assert_error, vec_of_owned};

  #[rustfmt::skip]
  #[rstest]
  #[case::default_command(         "--tree t.nwk --dates m.tsv",                  "timetree --tree t.nwk --metadata m.tsv")]
  #[case::default_after_root_flag( "-v --tree t.nwk",                             "-v timetree --tree t.nwk")]
  #[case::default_after_root_value("--verbosity info --tree t.nwk",               "--verbosity info timetree --tree t.nwk")]
  #[case::default_after_level(     "--verbose 3 --tree t.nwk",                    "--verbosity info timetree --tree t.nwk")]
  #[case::default_after_short_flag("-t t.nwk --dates m.tsv",                      "timetree -t t.nwk --metadata m.tsv")]
  #[case::explicit_command(        "timetree --tree t.nwk",                       "timetree --tree t.nwk")]
  #[case::no_arguments(            "",                                            "")]
  #[case::root_help(               "-h",                                          "-h")]
  #[case::root_version(            "--version",                                   "--version")]
  #[case::unknown_command(         "timetre --tree t.nwk",                        "timetre --tree t.nwk")]
  #[case::version_command(         "version",                                     "--version")]
  #[case::version_after_root_value("--color never version",                       "--color never --version")]
  #[trace]
  fn test_v0_args_command(#[case] argv: &str, #[case] expected: &str) {
    assert_eq!(rewritten_args(expected), rewrite(argv).unwrap().args);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::aln(                      "ancestral --aln a.fasta",                     "ancestral --alignment a.fasta")]
  #[case::aln_joined(               "ancestral --aln=a.fasta",                     "ancestral --alignment=a.fasta")]
  #[case::aln_in_prune(             "prune --aln a.fasta",                         "prune --alignment a.fasta")]
  #[case::aln_in_clock(             "clock --tree t.nwk --aln a.fasta",            "clock --tree t.nwk")]
  #[case::sequence_length_in_clock( "clock --sequence-length 1000 --tree t.nwk",   "clock --tree t.nwk")]
  #[case::sequence_length_timetree( "timetree --sequence-length 1000",             "timetree --sequence-length 1000")]
  #[case::dates(                    "clock --dates m.tsv",                         "clock --metadata m.tsv")]
  #[case::states(                   "mugration --states m.tsv",                    "mugration --metadata m.tsv")]
  #[case::name_column(              "timetree --name-column strain",               "timetree --metadata-id-columns strain")]
  #[case::outdir(                   "homoplasy --outdir out",                      "homoplasy --output-all out")]
  #[case::gtr(                      "ancestral --gtr JC69",                        "ancestral --model JC69")]
  #[case::gtr_params(               "timetree --gtr-params kappa=2 pis=0.1,0.2,0.3,0.4 --tree t.nwk",
                                                                                   "timetree --model-params kappa=2 pis=0.1,0.2,0.3,0.4 --tree t.nwk")]
  #[case::rng_seed(                 "homoplasy --rng-seed 7",                      "homoplasy --seed 7")]
  #[case::rng_seed_in_clock(        "clock --rng-seed 7",                          "clock")]
  #[case::rng_seed_in_mugration(    "mugration --rng-seed 7",                      "mugration")]
  #[case::aa(                       "ancestral --aa",                              "ancestral --alphabet aa")]
  #[case::aa_with_jtt92_gtr(        "ancestral --aa --gtr JTT92",                  "ancestral --alphabet aa-no-stop --model JTT92")]
  #[case::aa_with_jtt92_model(      "ancestral --model=jtt92 --aa",                "ancestral --model=jtt92 --alphabet aa-no-stop")]
  #[case::keep_overhangs(           "timetree --keep-overhangs",                   "timetree --gap-fill none")]
  #[case::marginal(                 "ancestral --marginal",                        "ancestral --method-anc marginal")]
  #[case::marginal_with_method(     "ancestral --marginal --method-anc fitch",     "ancestral --method-anc parsimony")]
  #[case::method_probabilistic(     "ancestral --method-anc probabilistic",        "ancestral --method-anc marginal")]
  #[case::method_ml(                "homoplasy --method-anc ml",                   "homoplasy --method-anc marginal")]
  #[case::method_fitch(             "homoplasy --method-anc fitch",                "homoplasy --method-anc parsimony")]
  #[case::method_parsimony(         "ancestral --method-anc parsimony",            "ancestral --method-anc parsimony")]
  #[case::timetree_method_ml(       "timetree --method-anc ml --tree t.nwk",       "timetree --tree t.nwk")]
  #[case::stochastic_resolve(       "timetree --stochastic-resolve",               "timetree")]
  #[case::reroot_best(              "timetree --reroot best",                      "timetree --reroot least-squares")]
  #[case::reroot_residual(          "clock --reroot residual",                     "clock --reroot least-squares")]
  #[case::reroot_res(               "optimize --reroot res",                       "optimize --reroot least-squares")]
  #[case::reroot_min_dev(           "timetree --reroot min_dev",                   "timetree --reroot min-dev")]
  #[case::reroot_oldest(            "timetree --reroot oldest",                    "timetree --reroot oldest")]
  #[case::reroot_least_squares(     "timetree --reroot=least-squares",             "timetree --reroot=least-squares")]
  #[case::reroot_tips(              "clock --reroot A B --tree t.nwk",             "clock --outgroup A B --tree t.nwk")]
  #[case::reroot_tip_joined(        "clock --reroot=A",                            "clock --outgroup=A")]
  #[case::reroot_tip_list_split_at_commas_which_unquoted_newick_names_lack("clock --reroot A,B", "clock --outgroup A,B")]
  #[case::coalescent_const(         "timetree --coalescent const",                 "timetree --coalescent opt")]
  #[case::coalescent_fixed(         "timetree --coalescent 0.01",                  "timetree --coalescent 0.01")]
  #[case::time_marginal_false(      "timetree --time-marginal false",              "timetree")]
  #[case::time_marginal_never(      "timetree --time-marginal=never",              "timetree")]
  #[case::verbose_0(                "clock --verbose 0",                           "clock --verbosity error")]
  #[case::verbose_1(                "clock --verbose 1",                           "clock --verbosity warn")]
  #[case::verbose_2(                "clock --verbose 2",                           "clock --verbosity info")]
  #[case::verbose_3(                "clock --verbose 3",                           "clock --verbosity info")]
  #[case::verbose_4(                "clock --verbose 4",                           "clock --verbosity debug")]
  #[case::verbose_6_joined(         "clock --verbose=6",                           "clock --verbosity debug")]
  #[case::verbose_count(            "clock --verbose --tree t.nwk",                "clock --verbose --tree t.nwk")]
  #[case::verbose_in_tool_command(  "schema --verbose 2 write",                    "schema --verbosity info write")]
  #[trace]
  fn test_v0_args_silent_rewrite(#[case] argv: &str, #[case] expected: &str) {
    assert_eq!(rewritten_args(expected), rewrite(argv).unwrap().args);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::after_separator(        "ancestral --tree t.nwk -- --aln a.fasta",     "ancestral --tree t.nwk -- --aln a.fasta")]
  #[case::joined_flag_value(      "ancestral --tree=--aln --alignment a.fasta",  "ancestral --tree=--aln --alignment a.fasta")]
  #[case::stdin_value(            "ancestral --tree - --aln a.fasta",            "ancestral --tree - --alignment a.fasta")]
  #[case::path_named_like_command("--dates version --tree arg",                  "timetree --metadata version --tree arg")]
  #[case::method_name_as_path(    "timetree --tree best --metadata min_dev",     "timetree --tree best --metadata min_dev")]
  #[case::zero_as_path(           "timetree --tree 0 --clock-filter 0",          "timetree --tree 0 --clock-filter 0")]
  #[case::short_attached_value(   "timetree -Oout -j4 --aln a.fasta",            "timetree -Oout -j4 --alignment a.fasta")]
  #[case::short_separate_value(   "timetree -O out --aln a.fasta",               "timetree -O out --alignment a.fasta")]
  #[case::level_after_values(     "timetree --alignment a.fasta b.fasta --verbose 2",
                                                                                 "timetree --alignment a.fasta b.fasta --verbosity info")]
  #[case::fixed_count(            "timetree --relax 1.0 0.5 --aln a.fasta",      "timetree --relax 1.0 0.5 --alignment a.fasta")]
  #[case::flag_of_other_command(  "mugration --aln a.fasta",                     "mugration --aln a.fasta")]
  #[trace]
  fn test_v0_args_reads_values_as_clap_does(#[case] argv: &str, #[case] expected: &str) {
    assert_eq!(rewritten_args(expected), rewrite(argv).unwrap().args);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(             "--coalescent 0",                         "")]
  #[case::zero_float(       "--coalescent 0.0",                       "")]
  #[case::zero_exponent(    "--coalescent=0e0",                       "")]
  #[case::negative_zero(    "--coalescent -0.0",                      "")]
  #[case::negative_zero_joined("--coalescent=-0.0",                   "")]
  #[case::last_zero_wins(   "--coalescent opt --coalescent 0",        "")]
  #[case::last_value_wins(  "--coalescent 0 --coalescent skyline",    "--coalescent skyline")]
  #[trace]
  fn test_v0_args_coalescent_zero(#[case] coalescent: &str, #[case] expected: &str) {
    let actual = rewrite(&format!("timetree --tree t.nwk {coalescent}")).unwrap();
    let cleared = if expected.is_empty() { btreeset! {"coalescent"} } else { btreeset! {} };
    assert_eq!(
      (rewritten_args(&format!("timetree --tree t.nwk {expected}")), cleared),
      (actual.args, actual.cleared)
    );
  }

  #[test]
  fn test_v0_args_returns_the_messages_of_the_silent_rewrites() {
    let actual = rewrite("--aln a.fasta --stochastic-resolve --coalescent 0 --verbose 2").unwrap();
    assert_eq!(
      V0Rewrite {
        args: rewritten_args("timetree --alignment a.fasta --verbosity info"),
        messages: vec_of_owned![
          "TreeTime 0.x command line without a command read as 'timetree'",
          "TreeTime 0.x argument '--aln a.fasta' read as '--alignment a.fasta'",
          "TreeTime 0.x argument '--stochastic-resolve' ignored: TreeTime 1 resolves polytomies stochastically",
          "TreeTime 0.x argument '--coalescent 0' clears the setting 'coalescent': TreeTime 0.x reads 0 as no coalescent prior",
          "TreeTime 0.x argument '--verbose 2' read as '--verbosity info'",
        ],
        cleared: btreeset! {"coalescent"},
      },
      actual
    );
  }

  #[test]
  fn test_v0_args_leaves_a_treetime_1_command_line_unchanged() {
    let argv =
      "timetree --tree t.nwk --metadata m.tsv --alignment a.fasta --reroot min-dev --coalescent skyline -O out";
    assert_eq!(
      V0Rewrite {
        args: rewritten_args(argv),
        messages: vec![],
        cleared: btreeset! {},
      },
      rewrite(argv).unwrap()
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reroot_ml(            "timetree --reroot ML",                   "--reroot=ML is not available in TreeTime 1: this method uses the covariation model of TreeTime 0.x, which TreeTime 1 removed; use --reroot least-squares")]
  #[case::reroot_min_dev_ml(    "clock --reroot min_dev_ML",              "--reroot=min_dev_ML is not available in TreeTime 1: this method uses the covariation model of TreeTime 0.x, which TreeTime 1 removed; use --reroot min-dev")]
  #[case::time_marginal_true(   "timetree --time-marginal true",          "--time-marginal=true is not available in TreeTime 1: use --confidence")]
  #[case::time_marginal_assign( "timetree --time-marginal assign",        "--time-marginal=assign is not available in TreeTime 1: use --confidence")]
  #[case::time_marginal_always( "--time-marginal always",                 "--time-marginal=always is not available in TreeTime 1: use --confidence")]
  #[case::time_marginal_final(  "timetree --time-marginal=only-final",    "--time-marginal=only-final is not available in TreeTime 1: use --confidence")]
  #[case::joint_branch_lengths( "timetree --branch-length-mode joint",    "--branch-length-mode=joint is not available in TreeTime 1: TreeTime 1 estimates branch lengths by marginal likelihood; use marginal or auto")]
  #[case::timetree_parsimony(   "timetree --method-anc parsimony",        "--method-anc=parsimony is not available in TreeTime 1: timetree reconstructs ancestral sequences by marginal likelihood only")]
  #[case::timetree_fitch(       "timetree --method-anc fitch",            "--method-anc=fitch is not available in TreeTime 1: timetree reconstructs ancestral sequences by marginal likelihood only")]
  #[case::greedy_resolve(       "timetree --greedy-resolve",              "--greedy-resolve is not available in TreeTime 1: TreeTime 1 resolves polytomies stochastically only")]
  #[case::covariation(          "clock --covariation",                    "--covariation is not available in TreeTime 1: TreeTime 1 fits the root-to-tip regression by ordinary least squares, without the covariation model")]
  #[case::tip_slack(            "timetree --tip-slack 3",                 "--tip-slack is not available in TreeTime 1: TreeTime 1 fits the root-to-tip regression by ordinary least squares, without the covariation model")]
  #[case::clock_filter_method(  "clock --clock-filter-method residual",   "--clock-filter-method is not available in TreeTime 1")]
  #[case::prune_outliers(       "clock --prune-outliers",                 "--prune-outliers is not available in TreeTime 1")]
  #[case::plot_tree(            "timetree --plot-tree tree.pdf",          "--plot-tree is not available in TreeTime 1")]
  #[case::plot_rtt(             "clock --plot-rtt rtt.pdf",               "--plot-rtt is not available in TreeTime 1")]
  #[case::tip_labels(           "timetree --tip-labels",                  "--tip-labels is not available in TreeTime 1")]
  #[case::no_tip_labels(        "timetree --no-tip-labels",               "--no-tip-labels is not available in TreeTime 1")]
  #[case::n_branches_posterior( "timetree --n-branches-posterior",        "--n-branches-posterior is not available in TreeTime 1")]
  #[case::vcf_reference(        "ancestral --vcf-reference ref.fasta",    "--vcf-reference is not available in TreeTime 1")]
  #[case::custom_gtr(           "homoplasy --custom-gtr gtr.txt",         "--custom-gtr is not available in TreeTime 1")]
  #[case::arg_command(          "arg --trees a.nwk b.nwk",                "the TreeTime 0.x command 'arg' is not available in TreeTime 1")]
  #[case::every_message(        "clock --covariation --plot-rtt rtt.pdf", "--covariation is not available in TreeTime 1: TreeTime 1 fits the root-to-tip regression by ordinary least squares, without the covariation model\n--plot-rtt is not available in TreeTime 1")]
  #[trace]
  fn test_v0_args_message(#[case] argv: &str, #[case] expected: &str) {
    assert_error!(rewrite(argv), expected);
  }

  #[test]
  fn test_v0_args_targets_exist_in_their_commands() {
    let root = treetime_command();
    let expected: Vec<String> = vec![];

    let actual = V0_SPELLINGS
      .iter()
      .flat_map(|entry| entry.commands.iter().map(move |name| (entry, *name)))
      .flat_map(|(entry, name)| {
        let command = find_command(&root, name);
        entry.rules.iter().filter_map(move |rule| {
          let target = match rule.then {
            Then::Rename(flag) => (flag, None),
            Then::Replace(tokens) => (tokens[0], tokens.get(1).copied()),
            _ => return None,
          };
          (!accepts(command, target)).then(|| format!("{name} {}: {target:?}", entry.flag))
        })
      })
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_v0_args_spellings_without_values_are_no_treetime_1_flags() {
    let root = treetime_command();
    let expected: Vec<String> = vec![];

    let actual = V0_SPELLINGS
      .iter()
      .filter(|entry| entry.rules.iter().all(|rule| matches!(rule.values, Values::Any)))
      .flat_map(|entry| entry.commands.iter().map(move |name| (entry.flag, *name)))
      .filter(|(flag, name)| long_flag(find_command(&root, name), flag).is_some())
      .map(|(flag, name)| format!("{name} {flag}"))
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_v0_args_list_each_spelling_once_per_command() {
    let expected: Vec<(&str, &str)> = vec![];

    let actual = V0_SPELLINGS
      .iter()
      .flat_map(|entry| entry.commands.iter().map(move |name| (*name, entry.flag)))
      .duplicates()
      .collect_vec();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_v0_args_every_command_lists_every_command() {
    let root = treetime_command();
    let expected: BTreeSet<&str> = iter::once(ROOT)
      .chain(
        root
          .get_subcommands()
          .map(Command::get_name)
          .filter(|name| *name != "help"),
      )
      .collect();

    let actual: BTreeSet<&str> = EVERY_COMMAND.iter().copied().collect();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_v0_args_reads_a_path_that_is_not_utf8_as_a_value() {
    let argv = vec![
      OsString::from("treetime"),
      OsString::from("ancestral"),
      OsString::from("--tree"),
      non_utf8(),
      OsString::from("--aln"),
      OsString::from("a.fasta"),
    ];
    let expected = vec![
      OsString::from("treetime"),
      OsString::from("ancestral"),
      OsString::from("--tree"),
      non_utf8(),
      OsString::from("--alignment"),
      OsString::from("a.fasta"),
    ];
    assert_eq!(expected, rewrite_v0_args(argv).unwrap().args);
  }

  #[test]
  fn test_v0_args_reads_an_attached_path_that_is_not_utf8_as_a_value() {
    let argv = vec![
      OsString::from("treetime"),
      OsString::from("ancestral"),
      OsString::from("--tree"),
      OsString::from("t.nwk"),
      attached("--aln=", &non_utf8()),
    ];
    let expected = vec![
      OsString::from("treetime"),
      OsString::from("ancestral"),
      OsString::from("--tree"),
      OsString::from("t.nwk"),
      attached("--alignment=", &non_utf8()),
    ];
    assert_eq!(expected, rewrite_v0_args(argv).unwrap().args);
  }

  mod helpers {
    use crate::cli::treetime_cli::TreetimeArgs;
    use crate::cli::v0_args::{ROOT, V0Rewrite, rewrite_v0_args};
    use clap::{Arg, Command, CommandFactory};
    use eyre::Report;
    use std::ffi::OsString;
    use std::iter;
    use std::os::unix::ffi::OsStringExt;

    pub(crate) fn rewrite(argv: &str) -> Result<V0Rewrite, Report> {
      rewrite_v0_args(rewritten_args(argv))
    }

    pub(crate) fn rewritten_args(argv: &str) -> Vec<OsString> {
      iter::once("treetime")
        .chain(argv.split_whitespace())
        .map(OsString::from)
        .collect()
    }

    pub(crate) fn treetime_command() -> Command {
      let mut command = TreetimeArgs::command();
      command.build();
      command
    }

    pub(crate) fn find_command<'a>(root: &'a Command, name: &str) -> &'a Command {
      if name == ROOT {
        root
      } else {
        root.find_subcommand(name).unwrap()
      }
    }

    pub(crate) fn long_flag<'a>(command: &'a Command, flag: &str) -> Option<&'a Arg> {
      let name = flag.strip_prefix("--").unwrap();
      command.get_arguments().find(|arg| arg.get_long() == Some(name))
    }

    pub(crate) fn accepts(command: &Command, (flag, value): (&str, Option<&str>)) -> bool {
      let Some(arg) = long_flag(command, flag) else {
        return false;
      };
      let possible = arg.get_possible_values();
      value.is_none_or(|value| possible.is_empty() || possible.iter().any(|known| known.matches(value, false)))
    }

    pub(crate) fn attached(flag: &str, value: &OsString) -> OsString {
      let mut token = OsString::from(flag);
      token.push(value);
      token
    }

    pub(crate) fn non_utf8() -> OsString {
      OsString::from_vec(vec![b't', 0xff, b'.', b'n', b'w', b'k'])
    }
  }
}
