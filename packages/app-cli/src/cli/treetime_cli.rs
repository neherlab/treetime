use crate::cli::color::ColorArgs;
use crate::cli::config::overlay_config;
use crate::cli::schema::SchemaTarget;
use crate::cli::verbosity::Verbosity;
use app_commands::command::RawConfig;
use app_commands::commands::ancestral::args::TreetimeAncestralArgsRaw;
use app_commands::commands::clock::args::TreetimeClockArgsRaw;
use app_commands::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
use app_commands::commands::mugration::args::TreetimeMugrationArgsRaw;
use app_commands::commands::optimize::args::TreetimeOptimizeArgsRaw;
use app_commands::commands::prune::args::TreetimePruneArgsRaw;
use app_commands::commands::shared::jobs::JobsArgs;
use app_commands::commands::timetree::args::TreetimeTimetreeArgsRaw;
use app_commands::version::LONG_VERSION;
use clap::{ArgMatches, ColorChoice, CommandFactory, FromArgMatches, Parser, Subcommand, ValueEnum, ValueHint};
use clap_complete::{Shell, generate};
use clap_complete_fig::Fig;
use deser::Serialize;
use eyre::Report;
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

pub(crate) fn treetime_parse_cli_args<I, T>(
  argv: I,
  color: ColorChoice,
) -> Result<(TreetimeArgs, ArgMatches), clap::Error>
where
  I: IntoIterator<Item = T>,
  T: Into<OsString> + Clone,
{
  let matches = TreetimeArgs::command().color(color).try_get_matches_from(argv)?;
  let args = TreetimeArgs::from_arg_matches(&matches)?;
  Ok((args, matches))
}

pub(crate) fn resolve_args(mut args: TreetimeArgs, matches: &ArgMatches) -> Result<TreetimeArgs, Report> {
  if let Some((_, sub_matches)) = matches.subcommand() {
    resolve_command_config(&mut args.command, sub_matches)?;
  }
  Ok(args)
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

fn resolve_command_config(command: &mut TreetimeCommands, matches: &ArgMatches) -> Result<(), Report> {
  match command {
    TreetimeCommands::Timetree(args) => resolve(args.as_mut(), matches),
    TreetimeCommands::Optimize(args) => resolve(args, matches),
    TreetimeCommands::Prune(args) => resolve(args, matches),
    TreetimeCommands::Ancestral(args) => resolve(args, matches),
    TreetimeCommands::Clock(args) => resolve(args, matches),
    TreetimeCommands::Homoplasy(args) => resolve(args, matches),
    TreetimeCommands::Mugration(args) => resolve(args, matches),
    _ => Ok(()),
  }
}

#[derive(Subcommand, Debug, Serialize)]
#[clap(verbatim_doc_comment)]
#[deser(rename_all = "kebab-case")]
pub(crate) enum TreetimeCommands {
  /// Estimates time trees from an initial tree topology, a set of date constraints (e.g. tip dates), and an alignment (optional).
  Timetree(Box<TreetimeTimetreeArgsRaw>),

  /// Reconstructs ancestral sequences and maps mutations to the tree. The `--output-*` flags select the outputs, such as reconstructed sequences, trees, and the inferred GTR model. With marginal reconstruction, each branch of the Nexus trees, and of the Newick trees with `--newick-annotations beast` or `nhx`, carries its mutations as a comment, like `A45G,G136T`, with 1-based positions.
  Ancestral(TreetimeAncestralArgsRaw),

  /// Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. It will reroot the tree to maximize the clock-like signal and recalculate branch length unless run with `--keep-root`.
  Clock(TreetimeClockArgsRaw),

  /// Reconstructs discrete ancestral states, for example geographic location, host, or similar. In addition to ancestral states, a GTR model of state transitions is inferred.
  Mugration(TreetimeMugrationArgsRaw),

  /// Reconstructs ancestral sequences and maps mutations to the tree. The tree is then scanned for homoplasies. An excess number of homoplasies might suggest contamination, recombination, culture adaptation or similar.
  Homoplasy(TreetimeHomoplasyArgsRaw),

  /// Optimizes the branch lengths and likelihood of a phylogenetic tree given aligned sequences.
  Optimize(TreetimeOptimizeArgsRaw),

  /// Prunes short branches and/or branches without mutations from a phylogenetic tree.
  Prune(TreetimePruneArgsRaw),

  /// Runs an ordered list of analysis commands from one config file, on one dataset, in one process.
  ///
  /// The config file (JSON or YAML) lists named steps, each an analysis command with its arguments.
  /// Steps run sequentially; a step may reference an earlier step's outputs with
  /// `{{ steps.<name>.outputs.<selection> }}`. Use `--dry-run` to print the resolved plan without
  /// running anything, and `--steps` to run only a subset.
  Pipeline(TreetimePipelineArgs),

  /// Write JSON Schema definitions for TreeTime data types
  Schema {
    #[clap(subcommand)]
    command: TreetimeSchemaCommands,
  },

  /// Generate shell completions.
  ///
  /// This will print the completions file contents to the console. Refer to your shell's documentation on how to install the completions.
  ///
  /// Example for Ubuntu Linux:
  ///
  ///    treetime completions bash > ~/.local/share/bash-completion/treetime
  ///
  Completions {
    /// Name of the shell to generate appropriate completions
    #[clap(value_name = "SHELL", default_value_t = String::from("bash"), value_parser = SHELLS.clone())]
    shell: String,
  },

  /// Print CLI reference documentation in Markdown format
  #[clap(name = "help-markdown")]
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
  #[clap(long, value_hint = ValueHint::FilePath)]
  pub config: PathBuf,

  /// Folder of the outputs: each step writes into a subfolder named after the step. Takes precedence over the
  /// top-level `output_all` of the pipeline file. A relative path resolves from the working directory.
  #[clap(long, short = 'O', value_hint = ValueHint::DirPath)]
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
}

#[derive(Subcommand, Debug, Serialize)]
#[deser(rename_all = "kebab-case")]
pub(crate) enum TreetimeSchemaCommands {
  /// Write the JSON schema of a config file or of all TreeTime data types.
  ///
  /// For example:
  ///
  ///    treetime schema write --for timetree -o timetree.schema.json
  Write(TreetimeSchemaArgs),
}

#[derive(Parser, Debug, Serialize)]
pub(crate) struct TreetimeSchemaArgs {
  /// Which schema to generate
  #[clap(long = "for", value_enum, default_value_t = SchemaTarget::default())]
  pub target: SchemaTarget,

  /// Output file or directory (use "-" for stdout). Directory required when --for=all
  #[clap(short = 'o', long, value_hint = ValueHint::AnyPath)]
  pub output: Option<PathBuf>,
}

fn resolve<T: RawConfig>(args: &mut T, matches: &ArgMatches) -> Result<(), Report> {
  overlay_config(args, matches)
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
