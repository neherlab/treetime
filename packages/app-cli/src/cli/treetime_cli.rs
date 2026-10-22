use crate::cli::color::ColorArgs;
use crate::cli::command_help::{
  ANCESTRAL_EXAMPLES, ANCESTRAL_USAGE, CLOCK_EXAMPLES, CLOCK_USAGE, COMPLETIONS_EXAMPLES, COMPLETIONS_USAGE, HELP_HINT,
  HELP_HINT_OUTPUTS, HELP_MARKDOWN_EXAMPLES, HELP_MARKDOWN_USAGE, HOMOPLASY_EXAMPLES, HOMOPLASY_USAGE, MUGRATION_EXAMPLES, MUGRATION_USAGE, OPTIMIZE_EXAMPLES, OPTIMIZE_USAGE,
  PIPELINE_EXAMPLES, PIPELINE_USAGE, PRUNE_EXAMPLES, PRUNE_USAGE, SCHEMA_WRITE_EXAMPLES, SCHEMA_WRITE_USAGE,
  TIMETREE_EXAMPLES, TIMETREE_USAGE, after_help, after_long_help,
};
use crate::cli::config::overlay_config;
use crate::cli::schema::SchemaTarget;
use crate::cli::v0_args::{V0Rewrite, rewrite_v0_args};
use crate::cli::verbosity::Verbosity;
use app_commands::command::RawConfig;
use app_commands::command_config::CommandConfig;
use app_commands::commands::ancestral::args::{ANCESTRAL_ABOUT, TreetimeAncestralArgsRaw};
use app_commands::commands::clock::args::{CLOCK_ABOUT, TreetimeClockArgsRaw};
use app_commands::commands::homoplasy::args::{HOMOPLASY_ABOUT, TreetimeHomoplasyArgsRaw};
use app_commands::commands::mugration::args::{MUGRATION_ABOUT, TreetimeMugrationArgsRaw};
use app_commands::commands::optimize::args::{OPTIMIZE_ABOUT, TreetimeOptimizeArgsRaw};
use app_commands::commands::prune::args::{PRUNE_ABOUT, TreetimePruneArgsRaw};
use app_commands::commands::shared::jobs::JobsArgs;
use app_commands::commands::timetree::args::{TIMETREE_ABOUT, TreetimeTimetreeArgsRaw};
use app_commands::version::LONG_VERSION;
use clap::{ArgMatches, ColorChoice, CommandFactory, FromArgMatches, Parser, Subcommand, ValueEnum, ValueHint};
use clap_complete::{Shell, generate};
use clap_complete_fig::Fig;
use deser::Serialize;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeSet;
use std::ffi::OsString;
use std::fmt::Debug;
use std::io;
use std::path::PathBuf;
use std::sync::LazyLock;
use treetime_utils::init::clap_styles::styles;
use treetime_utils::make_report;

static SHELLS: LazyLock<Vec<&'static str>> =
  LazyLock::new(|| vec!["bash", "elvish", "fish", "fig", "powershell", "zsh"]);

pub(crate) fn generate_shell_completions(shell: &str) -> Result<(), Report> {
  let mut command = TreetimeArgs::command();

  if shell.to_lowercase() == "fig" {
    generate(Fig, &mut command, "treetime", &mut io::stdout());
    return Ok(());
  }

  let generator = <Shell as ValueEnum>::from_str(&shell.to_lowercase(), true)
    .map_err(|err| make_report!("{err}: Possible values: {}", SHELLS.join(", ")))?;

  let bin_name = command.get_name().to_owned();

  generate(generator, &mut command, bin_name, &mut io::stdout());

  Ok(())
}

pub(crate) fn treetime_parse_cli_args<I, T>(argv: I, color: ColorChoice) -> Result<ParsedCli, CliParseError>
where
  I: IntoIterator<Item = T>,
  T: Into<OsString> + Clone,
{
  let V0Rewrite {
    args: argv,
    messages: v0_messages,
    cleared,
  } = rewrite_v0_args(argv.into_iter().map(Into::into).collect_vec()).map_err(CliParseError::Argument)?;
  let matches = TreetimeArgs::command()
    .color(color)
    .try_get_matches_from(argv)
    .map_err(CliParseError::Clap)?;
  let args = TreetimeArgs::from_arg_matches(&matches).map_err(CliParseError::Clap)?;
  Ok(ParsedCli {
    args,
    matches,
    v0_messages,
    cleared,
  })
}

pub(crate) fn resolve_args(parsed: ParsedCli) -> Result<TreetimeArgs, Report> {
  let ParsedCli {
    mut args,
    matches,
    cleared,
    ..
  } = parsed;
  if let Some((_, sub_matches)) = matches.subcommand() {
    resolve_command_config(&mut args.command, sub_matches, &cleared)?;
  }
  Ok(args)
}

pub(crate) struct ParsedCli {
  pub args: TreetimeArgs,
  pub matches: ArgMatches,
  pub v0_messages: Vec<String>,
  pub cleared: BTreeSet<&'static str>,
}

#[derive(Debug)]
pub(crate) enum CliParseError {
  Clap(clap::Error),
  Argument(Report),
}

#[derive(Parser, Debug, Serialize)]
#[clap(name = "treetime")]
#[clap(author, version = LONG_VERSION)]
#[clap(verbatim_doc_comment)]
#[clap(styles = styles())]
#[clap(after_help = "Exit status: 0 on success, 1 on failure, 2 on a usage error (invalid flags or config file)")]
/// Maximum-likelihood phylodynamic inference
///
/// Repository:  https://github.com/neherlab/treetime
/// Publication: https://academic.oup.com/ve/article/4/1/vex042/4794731
pub(crate) struct TreetimeArgs {
  #[clap(flatten, next_help_heading = "Verbosity")]
  pub verbosity: Verbosity,

  #[clap(flatten, next_help_heading = "Verbosity")]
  pub color: ColorArgs,

  #[clap(subcommand)]
  pub command: TreetimeCommands,
}

fn resolve_command_config(
  command: &mut TreetimeCommands,
  matches: &ArgMatches,
  cleared: &BTreeSet<&str>,
) -> Result<(), Report> {
  match command {
    TreetimeCommands::Timetree(args) => overlay_config(args.as_mut(), matches, cleared),
    TreetimeCommands::Optimize(args) => overlay_config(args, matches, cleared),
    TreetimeCommands::Prune(args) => overlay_config(args, matches, cleared),
    TreetimeCommands::Ancestral(args) => overlay_config(args, matches, cleared),
    TreetimeCommands::Clock(args) => overlay_config(args, matches, cleared),
    TreetimeCommands::Homoplasy(args) => overlay_config(args, matches, cleared),
    TreetimeCommands::Mugration(args) => overlay_config(args, matches, cleared),
    _ => Ok(()),
  }
}

#[derive(Subcommand, Debug, Serialize)]
#[clap(verbatim_doc_comment)]
#[deser(rename_all = "kebab-case")]
pub(crate) enum TreetimeCommands {
  #[command(
    about = TIMETREE_ABOUT,
    override_usage = TIMETREE_USAGE,
    after_help = after_help(TIMETREE_EXAMPLES, HELP_HINT_OUTPUTS),
    after_long_help = after_long_help(TIMETREE_EXAMPLES)
  )]
  Timetree(Box<TreetimeTimetreeArgsRaw>),

  #[command(
    about = ANCESTRAL_ABOUT,
    long_about = format!(
      "{ANCESTRAL_ABOUT}\n\nWith marginal reconstruction, each branch of the Nexus trees, and of the Newick trees with \
       `--newick-annotations beast` or `nhx`, carries its mutations as a comment, like `A45G,G136T`, with 1-based \
       positions."
    ),
    override_usage = ANCESTRAL_USAGE,
    after_help = after_help(ANCESTRAL_EXAMPLES, HELP_HINT_OUTPUTS),
    after_long_help = after_long_help(ANCESTRAL_EXAMPLES)
  )]
  Ancestral(TreetimeAncestralArgsRaw),

  #[command(
    about = CLOCK_ABOUT,
    long_about = format!(
      "{CLOCK_ABOUT}\n\nThe regression quantifies how clock-like the tree is. Unless `--keep-root` is given, the run \
       reroots the tree to maximize the clock signal."
    ),
    override_usage = CLOCK_USAGE,
    after_help = after_help(CLOCK_EXAMPLES, HELP_HINT_OUTPUTS),
    after_long_help = after_long_help(CLOCK_EXAMPLES)
  )]
  Clock(TreetimeClockArgsRaw),

  #[command(
    about = MUGRATION_ABOUT,
    override_usage = MUGRATION_USAGE,
    after_help = after_help(MUGRATION_EXAMPLES, HELP_HINT_OUTPUTS),
    after_long_help = after_long_help(MUGRATION_EXAMPLES)
  )]
  Mugration(TreetimeMugrationArgsRaw),

  #[command(
    about = HOMOPLASY_ABOUT,
    long_about = format!(
      "{HOMOPLASY_ABOUT}\n\nAn excess of homoplasies can point to contamination, recombination, culture adaptation, or \
       similar processes."
    ),
    override_usage = HOMOPLASY_USAGE,
    after_help = after_help(HOMOPLASY_EXAMPLES, HELP_HINT_OUTPUTS),
    after_long_help = after_long_help(HOMOPLASY_EXAMPLES)
  )]
  Homoplasy(TreetimeHomoplasyArgsRaw),

  #[command(
    about = OPTIMIZE_ABOUT,
    override_usage = OPTIMIZE_USAGE,
    after_help = after_help(OPTIMIZE_EXAMPLES, HELP_HINT_OUTPUTS),
    after_long_help = after_long_help(OPTIMIZE_EXAMPLES)
  )]
  Optimize(TreetimeOptimizeArgsRaw),

  #[command(
    about = PRUNE_ABOUT,
    override_usage = PRUNE_USAGE,
    after_help = after_help(PRUNE_EXAMPLES, HELP_HINT_OUTPUTS),
    after_long_help = after_long_help(PRUNE_EXAMPLES)
  )]
  Prune(TreetimePruneArgsRaw),

  /// Run several analysis commands in order from one config file
  ///
  /// The config file (JSON or YAML) lists named steps, each an analysis command with its arguments. The steps run one
  /// after another on one dataset, in one process; a step may reference the outputs of an earlier step with
  /// `{{ steps.<name>.outputs.<selection> }}`. Use `--dry-run` to print the resolved plan without running anything,
  /// and `--steps` to run only some steps.
  #[command(
    override_usage = PIPELINE_USAGE,
    after_help = after_help(PIPELINE_EXAMPLES, HELP_HINT),
    after_long_help = after_long_help(PIPELINE_EXAMPLES)
  )]
  Pipeline(TreetimePipelineArgs),

  /// Write the JSON schemas of the config files and of the TreeTime data types
  Schema {
    #[clap(subcommand)]
    command: TreetimeSchemaCommands,
  },

  /// Print the shell completions of TreeTime
  ///
  /// Write the output to the completions file of your shell; the documentation of the shell names its location.
  #[command(
    override_usage = COMPLETIONS_USAGE,
    after_help = after_help(COMPLETIONS_EXAMPLES, HELP_HINT),
    after_long_help = after_long_help(COMPLETIONS_EXAMPLES)
  )]
  Completions {
    /// Name of the shell to generate appropriate completions
    #[clap(value_name = "SHELL", default_value_t = String::from("bash"), value_parser = SHELLS.clone())]
    shell: String,
  },

  /// Print the command-line reference in Markdown
  #[clap(name = "help-markdown")]
  #[command(
    override_usage = HELP_MARKDOWN_USAGE,
    after_help = after_help(HELP_MARKDOWN_EXAMPLES, HELP_HINT),
    after_long_help = after_long_help(HELP_MARKDOWN_EXAMPLES)
  )]
  HelpMarkdown,

  /// Print system information for debugging
  #[clap(hide = true)]
  Debug,
}

#[derive(Parser, Debug, Serialize)]
pub(crate) struct TreetimePipelineArgs {
  /// Pipeline configuration file (JSON or YAML) describing an ordered list of steps.
  ///
  /// A relative path in the file resolves from the folder of the file, and from the working directory when the file
  /// is read from standard input (`-`).
  #[clap(long, value_hint = ValueHint::FilePath, value_name = "FILE")]
  pub config: PathBuf,

  /// Folder of the outputs: each step writes into a subfolder named after the step. Takes precedence over the
  /// top-level `output_all` of the pipeline file. A relative path resolves from the working directory.
  #[clap(long, short = 'O', value_hint = ValueHint::DirPath, value_name = "DIR")]
  pub output_all: Option<PathBuf>,

  /// Run only these named steps, in the pipeline's list order.
  ///
  /// A referenced upstream step that is not selected must already have its outputs on disk.
  #[clap(long, num_args = 1.., value_delimiter = ',')]
  pub steps: Vec<String>,

  /// Resolve and print the plan (steps, inputs, outputs) without running any step.
  #[clap(long)]
  pub dry_run: bool,

  #[clap(flatten)]
  pub jobs_args: JobsArgs,
}

impl TreetimeCommands {
  pub(crate) fn jobs(&self) -> Option<usize> {
    match self {
      Self::Timetree(args) => Some(args.jobs_args.jobs),
      Self::Optimize(args) => Some(args.jobs_args.jobs),
      Self::Prune(args) => Some(args.jobs_args.jobs),
      Self::Ancestral(args) => Some(args.jobs_args.jobs),
      Self::Clock(args) => Some(args.jobs_args.jobs),
      Self::Homoplasy(args) => Some(args.jobs_args.jobs),
      Self::Mugration(args) => Some(args.jobs_args.jobs),
      Self::Pipeline(args) => Some(args.jobs_args.jobs),
      Self::Completions { .. } | Self::Schema { .. } | Self::HelpMarkdown | Self::Debug => None,
    }
  }

  pub(crate) fn dry_run_config(&self) -> Option<CommandConfig> {
    match self {
      Self::Timetree(args) => dry_run_config_of(args.as_ref(), args.config_args.dry_run),
      Self::Optimize(args) => dry_run_config_of(args, args.config_args.dry_run),
      Self::Prune(args) => dry_run_config_of(args, args.config_args.dry_run),
      Self::Ancestral(args) => dry_run_config_of(args, args.config_args.dry_run),
      Self::Clock(args) => dry_run_config_of(args, args.config_args.dry_run),
      Self::Homoplasy(args) => dry_run_config_of(args, args.config_args.dry_run),
      Self::Mugration(args) => dry_run_config_of(args, args.config_args.dry_run),
      Self::Pipeline(_)
      | Self::Completions { .. }
      | Self::Schema { .. }
      | Self::HelpMarkdown
      | Self::Debug => None,
    }
  }
}

fn dry_run_config_of<T: RawConfig>(args: &T, dry_run: bool) -> Option<CommandConfig> {
  dry_run.then(|| args.clone().into_config())
}

#[derive(Subcommand, Debug, Serialize)]
#[deser(rename_all = "kebab-case")]
pub(crate) enum TreetimeSchemaCommands {
  /// Write the JSON schema of a config file, or of every config file and TreeTime data type
  #[command(
    override_usage = SCHEMA_WRITE_USAGE,
    after_help = after_help(SCHEMA_WRITE_EXAMPLES, HELP_HINT),
    after_long_help = after_long_help(SCHEMA_WRITE_EXAMPLES)
  )]
  Write(TreetimeSchemaArgs),
}

#[derive(Parser, Debug, Serialize)]
pub(crate) struct TreetimeSchemaArgs {
  /// Config file or data type of the schema; `all` writes every schema into the folder of `-o`
  #[clap(long = "for", value_enum, value_name = "TARGET", default_value_t = SchemaTarget::default())]
  pub target: SchemaTarget,

  /// Output file, or `-` for standard output; with `--for all`, the output folder (default: the working folder)
  #[clap(short = 'o', long, value_hint = ValueHint::AnyPath, value_name = "PATH")]
  pub output: Option<PathBuf>,

  /// Print the files that the command would write, then stop without writing them
  #[clap(long)]
  pub dry_run: bool,
}

#[cfg(test)]
mod tests {
  use super::TreetimeArgs;
  use app_commands::commands::timetree::args::TreetimeTimetreeArgsRaw;
  use clap::error::ErrorKind;
  use clap::{Command, CommandFactory, Parser};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_treetime_cli_root_lists_the_analysis_commands_first() {
    let command = TreetimeArgs::command();
    let actual: Vec<&str> = command
      .get_subcommands()
      .filter(|subcommand| !subcommand.is_hide_set())
      .map(Command::get_name)
      .collect();
    assert_eq!(
      vec![
        "timetree",
        "ancestral",
        "clock",
        "mugration",
        "homoplasy",
        "optimize",
        "prune",
        "pipeline",
        "schema",
        "completions",
        "help-markdown",
      ],
      actual
    );
  }

  fn parse_timetree(extra: &[&str]) -> Result<TreetimeTimetreeArgsRaw, clap::Error> {
    let argv = std::iter::once("timetree").chain(extra.iter().copied());
    TreetimeTimetreeArgsRaw::try_parse_from(argv)
  }

  #[test]
  fn test_treetime_cli_coalescent_settings_are_unset_and_resolve_to_their_documented_defaults()
  -> Result<(), clap::Error> {
    let args = parse_timetree(&[])?;
    assert_eq!(
      (None, None, None, None),
      (
        args.n_skyline,
        args.skyline_stiffness,
        args.coalescent_confidence,
        args.gen_per_year
      )
    );
    assert_eq!(20, args.n_skyline());
    pretty_assert_ulps_eq!(2.0, args.skyline_stiffness());
    pretty_assert_ulps_eq!(2.0, args.coalescent_confidence());
    pretty_assert_ulps_eq!(50.0, args.gen_per_year());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::confidence(  &["--coalescent-confidence", "3.5"], (Some(3.5), None,    None))]
  #[case::n_skyline(   &["--n-skyline",             "4"],   (None,      Some(4), None))]
  #[case::gen_per_year(&["--gen-per-year",          "30"],  (None,      None,    Some(30.0)))]
  #[trace]
  fn test_treetime_cli_coalescent_args_parse(
    #[case] extra: &[&str],
    #[case] expected: (Option<f64>, Option<usize>, Option<f64>),
  ) {
    let args = parse_timetree(extra).unwrap();
    assert_eq!(expected, (args.coalescent_confidence, args.n_skyline, args.gen_per_year));
  }

  #[test]
  fn test_treetime_cli_rejects_the_skyline_n_points_flag() {
    match parse_timetree(&["--skyline-n-points", "5"]) {
      Ok(_) => panic!("--skyline-n-points must be rejected"),
      Err(err) => assert_eq!(ErrorKind::UnknownArgument, err.kind()),
    }
  }
}
