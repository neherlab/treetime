#![allow(unused_qualifications)]

use crate::cli::jobs::Jobs;
use crate::cli::verbosity::Verbosity;
use clap::{CommandFactory, Parser, Subcommand, ValueEnum};
use clap_complete::{Shell, generate};
use clap_complete_fig::Fig;
use eyre::{Report, eyre};
use lazy_static::lazy_static;
use std::fmt::Debug;
use std::io;
use treetime::commands::ancestral::anc_args::TreetimeAncestralArgs;
use treetime::commands::clock::clock_args::TreetimeClockArgs;
use treetime::commands::homoplasy::homoplasy_args::TreetimeHomoplasyArgs;
use treetime::commands::mugration::mugration_args::TreetimeMugrationArgs;
use treetime::commands::optimize::args::TreetimeOptimizeArgs;
use treetime::commands::timetree::timetree_args::TreetimeTimetreeArgs;
use treetime::utils::clap_styles::styles;
use treetime::utils::global_init::setup_logger;

lazy_static! {
  pub static ref SHELLS: Vec<&'static str> = ["bash", "elvish", "fish", "fig", "powershell", "zsh"].to_vec();
}

#[derive(Parser, Debug)]
#[clap(name = "treetime")]
#[clap(author, version)]
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

#[derive(Subcommand, Debug)]
#[clap(verbatim_doc_comment)]
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
  Timetree(TreetimeTimetreeArgs),

  // TODO: explain what this command does
  Optimize(TreetimeOptimizeArgs),

  /// Reconstructs ancestral sequences and maps mutations to the tree. The output consists of a file 'ancestral.fasta' with ancestral sequences and a tree 'annotated_tree.nexus' with mutations added as comments like A45G,G136T,..., number in SNPs used 1-based index by default. The inferred GTR model is written to stdout.
  Ancestral(TreetimeAncestralArgs),

  /// Calculates the root-to-tip regression and quantifies the 'clock-i-ness' of the tree. It will reroot the tree to maximize the clock-like signal and recalculate branch length unless run with --keep_root.
  Clock(TreetimeClockArgs),

  /// Reconstructs ancestral sequences and maps mutations to the tree. The tree is then scanned for homoplasies. An excess number of homoplasies might suggest contamination, recombination, culture adaptation or similar.
  Homoplasy(TreetimeHomoplasyArgs),

  /// Reconstructs discrete ancestral states, for example geographic location, host, or similar. In addition to ancestral states, a GTR model of state transitions is inferred.
  Mugration(TreetimeMugrationArgs),

  /// Estimates ancestral reassortment graph (ARG).
  Arg(TreetimeAncestralReassortmentGraphArgs),

  /// Print system information for debugging
  #[clap(hide = true)]
  Debug,
}

#[derive(Parser, Debug)]
pub struct TreetimeAncestralReassortmentGraphArgs;

pub fn generate_shell_completions(shell: &str) -> Result<(), Report> {
  let mut command = TreetimeArgs::command();

  if shell.to_lowercase() == "fig" {
    generate(Fig, &mut command, "treetime", &mut io::stdout());
    return Ok(());
  }

  let generator = <Shell as ValueEnum>::from_str(&shell.to_lowercase(), true)
    .map_err(|err| eyre!("{}: Possible values: {}", err, SHELLS.join(", ")))?;

  let bin_name = command.get_name().to_owned();

  generate(generator, &mut command, bin_name, &mut io::stdout());

  Ok(())
}

pub fn treetime_parse_cli_args() -> Result<TreetimeArgs, Report> {
  let args = TreetimeArgs::parse();
  setup_logger(args.verbosity.get_filter_level());
  Ok(args)
}
