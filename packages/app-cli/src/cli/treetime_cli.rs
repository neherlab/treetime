use crate::cli::config::overlay_config;
use crate::cli::jobs::Jobs;
use crate::cli::validate::Validate;
use crate::cli::verbosity::Verbosity;
use clap::{ArgMatches, CommandFactory, FromArgMatches, Parser, Subcommand, ValueEnum, ValueHint};
use clap_complete::{Shell, generate};
use clap_complete_fig::Fig;
use eyre::Report;
use serde::Serialize;
use serde::de::DeserializeOwned;
use std::fmt::Debug;
use std::io;
use std::path::PathBuf;
use std::sync::LazyLock;
use treetime::commands::ancestral::args::TreetimeAncestralArgs;
use treetime::commands::clock::args::TreetimeClockArgs;
use treetime::commands::homoplasy::args::TreetimeHomoplasyArgs;
use treetime::commands::mugration::args::TreetimeMugrationArgs;
use treetime::commands::optimize::args::TreetimeOptimizeArgs;
use treetime::commands::prune::args::TreetimePruneArgs;
use treetime::commands::timetree::args::TreetimeTimetreeArgs;
use treetime_schema::TreetimeSchemaFormat;
use treetime_utils::init::clap_styles::styles;
use treetime_utils::init::global::setup_logger;
use treetime_utils::make_report;

pub static SHELLS: LazyLock<Vec<&'static str>> =
  LazyLock::new(|| vec!["bash", "elvish", "fish", "fig", "powershell", "zsh"]);

#[derive(Parser, Debug, Serialize)]
#[clap(name = "treetime")]
#[clap(author, version = env!("TREETIME_LONG_VERSION"))]
#[clap(verbatim_doc_comment)]
#[clap(styles = styles())]
/// Maximum-likelihood phylodynamic inference
///
/// Documentation: https://treetime.readthedocs.io/en/stable/
/// Publication:   https://academic.oup.com/ve/article/4/1/vex042/4794731
pub struct TreetimeArgs {
  #[clap(flatten, next_help_heading = "Parallelism")]
  pub jobs: Jobs,

  #[clap(flatten, next_help_heading = "Verbosity")]
  pub verbosity: Verbosity,

  #[clap(subcommand)]
  pub command: TreetimeCommands,
}

#[derive(Subcommand, Debug, Serialize)]
#[clap(verbatim_doc_comment)]
#[serde(rename_all = "kebab-case")]
pub enum TreetimeCommands {
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

  /// Estimates time trees from an initial tree topology, a set of date constraints (e.g. tip dates), and an alignment (optional).
  Timetree(Box<TreetimeTimetreeArgs>),

  /// Optimizes the branch lengths and likelihood of a phylogenetic tree given aligned sequences.
  Optimize(TreetimeOptimizeArgs),

  /// Prunes short branches and/or branches without mutations from a phylogenetic tree.
  Prune(TreetimePruneArgs),

  /// Reconstructs ancestral sequences and maps mutations to the tree. The output consists of a file 'ancestral.fasta' with ancestral sequences and a tree 'ancestral.nexus' with mutations added as comments like A45G,G136T,..., number in SNPs used 1-based index by default. The inferred GTR model is written to stdout.
  Ancestral(TreetimeAncestralArgs),

  /// Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. It will reroot the tree to maximize the clock-like signal and recalculate branch length unless run with --keep_root.
  Clock(TreetimeClockArgs),

  /// Reconstructs ancestral sequences and maps mutations to the tree. The tree is then scanned for homoplasies. An excess number of homoplasies might suggest contamination, recombination, culture adaptation or similar.
  Homoplasy(TreetimeHomoplasyArgs),

  /// Reconstructs discrete ancestral states, for example geographic location, host, or similar. In addition to ancestral states, a GTR model of state transitions is inferred.
  Mugration(TreetimeMugrationArgs),

  /// Runs an ordered list of analysis commands from one config file, on one dataset, in one process.
  ///
  /// The config file (JSON or YAML) lists named steps, each an analysis command with its arguments.
  /// Steps run sequentially; a step may reference an earlier step's outputs with
  /// `{{ steps.<name>.outputs.<selection> }}`. Use `--check` to print the resolved plan without
  /// running anything, and `--steps` to run only a subset.
  Pipeline(TreetimePipelineArgs),

  /// Estimates ancestral reassortment graph (ARG).
  Arg(TreetimeAncestralReassortmentGraphArgs),

  /// Write JSON Schema definitions for TreeTime data types
  Schema(TreetimeSchemaArgs),

  /// Print CLI reference documentation in Markdown format
  #[clap(name = "help-markdown")]
  HelpMarkdown,

  /// Print system information for debugging
  #[clap(hide = true)]
  Debug,
}

#[derive(Parser, Debug, Serialize)]
pub struct TreetimePipelineArgs {
  /// Pipeline configuration file (JSON or YAML) describing an ordered list of steps.
  #[clap(long, value_hint = ValueHint::FilePath)]
  pub config: PathBuf,

  /// Run only these named steps (comma-separated), in the pipeline's list order.
  ///
  /// A referenced upstream step that is not selected must already have its outputs on disk.
  #[clap(long, value_delimiter = ',')]
  pub steps: Vec<String>,

  /// Resolve and print the plan (steps, inputs, outputs) without running any step.
  #[clap(long)]
  pub check: bool,
}

#[derive(Parser, Debug, Serialize)]
pub struct TreetimeSchemaArgs {
  /// Which schema to generate
  #[clap(long = "for", value_enum, default_value_t = TreetimeSchemaFormat::default())]
  pub for_format: TreetimeSchemaFormat,

  /// Output file or directory (use "-" for stdout). Directory required when --for=all
  #[clap(short = 'o', long, value_hint = ValueHint::AnyPath)]
  pub output: Option<PathBuf>,
}

#[derive(Parser, Debug, Serialize)]
pub struct TreetimeAncestralReassortmentGraphArgs;

pub fn generate_shell_completions(shell: &str) -> Result<(), Report> {
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

pub fn treetime_parse_cli_args() -> Result<TreetimeArgs, Report> {
  // `get_matches` (not `parse`) so we can consult per-argument value sources when merging `--config`.
  let matches = TreetimeArgs::command().get_matches();
  let mut args = TreetimeArgs::from_arg_matches(&matches)?;
  setup_logger(args.verbosity.get_filter_level());

  if let Some((_, sub_matches)) = matches.subcommand() {
    resolve_command_config(&mut args.command, sub_matches)?;
  }

  Ok(args)
}

/// Merge a `--config` file into the selected command's args, then validate the merged result.
///
/// Merge (`overlay_config`) and validation (`Validate`) are applied together so a value supplied only
/// by the config file both wins with the right precedence and satisfies required-argument checks.
/// Commands without a configuration object (completions, schema, help) have no `--config` flag and
/// are left untouched.
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

/// Overlay the `--config` file onto `args`, then validate the merged result.
fn resolve<T>(args: &mut T, matches: &ArgMatches) -> Result<(), Report>
where
  T: Serialize + DeserializeOwned + Default + Validate,
{
  overlay_config(args, matches)?;
  args.validate()
}

#[cfg(test)]
mod tests {
  use clap::Parser;
  use clap::error::ErrorKind;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::commands::timetree::args::TreetimeTimetreeArgs;
  use treetime_utils::pretty_assert_ulps_eq;

  // Timetree declares no clap-required arguments, so a bare invocation exercises defaults.
  fn parse_timetree(extra: &[&str]) -> Result<TreetimeTimetreeArgs, clap::Error> {
    let argv = std::iter::once("timetree").chain(extra.iter().copied());
    TreetimeTimetreeArgs::try_parse_from(argv)
  }

  #[test]
  fn test_treetime_cli_coalescent_defaults() -> Result<(), clap::Error> {
    // Defaults live in `TreetimeTimetreeArgs` (SmartDefault) and clap reads them via `default_value_t`.
    let args = parse_timetree(&[])?;
    pretty_assert_ulps_eq!(2.0, args.coalescent_confidence);
    assert_eq!(20, args.skyline_n_points);
    pretty_assert_ulps_eq!(50.0, args.gen_per_year);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::confidence(  &["--coalescent-confidence", "3.5"], 3.5, 20, 50.0)]
  #[case::n_points(    &["--skyline-n-points",      "4"],   2.0,  4, 50.0)]
  #[case::gen_per_year(&["--gen-per-year",          "30"],  2.0, 20, 30.0)]
  #[trace]
  fn test_treetime_cli_coalescent_args_parse(
    #[case] extra: &[&str],
    #[case] expected_confidence: f64,
    #[case] expected_n_points: usize,
    #[case] expected_gen_per_year: f64,
  ) {
    let args = parse_timetree(extra).unwrap();
    pretty_assert_ulps_eq!(expected_confidence, args.coalescent_confidence);
    assert_eq!(expected_n_points, args.skyline_n_points);
    pretty_assert_ulps_eq!(expected_gen_per_year, args.gen_per_year);
  }

  #[test]
  fn test_treetime_cli_skyline_n_points_rejects_below_two() {
    // The renamed `--skyline-n-points` keeps the >= 2 lower bound of the old `--n-skyline`.
    match parse_timetree(&["--skyline-n-points", "1"]) {
      Ok(_) => panic!("--skyline-n-points=1 must be rejected"),
      Err(err) => assert_eq!(ErrorKind::ValueValidation, err.kind()),
    }
  }

  #[test]
  fn test_treetime_cli_rejects_renamed_n_skyline_flag() {
    // The old `--n-skyline` spelling is gone; only `--skyline-n-points` is accepted.
    match parse_timetree(&["--n-skyline", "5"]) {
      Ok(_) => panic!("removed --n-skyline flag must be rejected"),
      Err(err) => assert_eq!(ErrorKind::UnknownArgument, err.kind()),
    }
  }
}
