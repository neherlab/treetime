#![allow(clippy::large_enum_variant)]
#![allow(clippy::struct_excessive_bools)]

use clap::{AppSettings, ArgEnum, CommandFactory, Parser, Subcommand};
use clap_complete::{generate, Generator, Shell};
use clap_complete_fig::Fig;
use clap_verbosity_flag::{Verbosity, WarnLevel};
use eyre::{eyre, Report};
use itertools::Itertools;
use lazy_static::lazy_static;
use log::LevelFilter;
use num_cpus;
use std::fmt::Debug;
use std::io;
use std::str::FromStr;
use treetime::ancestral::anc_args::TreetimeAncestralArgs;
use treetime::clock::clock_args::TreetimeClockArgs;
use treetime::homoplasy::homoplasy_args::TreetimeHomoplasyArgs;
use treetime::mugration::mugration_args::TreetimeMugrationArgs;
use treetime::timetree::timetree_args::TreetimeTimetreeArgs;
use treetime::utils::global_init::setup_logger;

lazy_static! {
  static ref SHELLS: &'static [&'static str] = &["bash", "elvish", "fish", "fig", "powershell", "zsh"];
  static ref VERBOSITIES: &'static [&'static str] = &["off", "error", "warn", "info", "debug", "trace"];
}

#[derive(Parser, Debug)]
#[clap(name = "treetime", trailing_var_arg = true)]
#[clap(author, version)]
#[clap(global_setting(AppSettings::DeriveDisplayOrder))]
#[clap(verbatim_doc_comment)]
/// Maximum-likelihood phylodynamic inference
///
/// Documentation: https://treetime.readthedocs.io/en/stable/
/// Publication:   https://academic.oup.com/ve/article/4/1/vex042/4794731
pub struct TreetimeArgs {
  #[clap(subcommand)]
  pub command: TreetimeCommands,

  /// Make output more quiet or more verbose
  #[clap(flatten)]
  pub verbose: Verbosity<WarnLevel>,

  /// Set verbosity level
  #[clap(long, global = true, conflicts_with = "verbose", conflicts_with = "silent", possible_values(VERBOSITIES.iter()))]
  pub verbosity: Option<LevelFilter>,

  /// Disable all console output. Same as --verbosity=off
  #[clap(long, global = true, conflicts_with = "verbose", conflicts_with = "verbosity")]
  pub silent: bool,

  /// Number of processing jobs. If not specified, all available CPU threads will be used.
  #[clap(global = true, long, short = 'j', default_value_t = num_cpus::get())]
  pub jobs: usize,
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
    #[clap(value_name = "SHELL", default_value_t = String::from("bash"), possible_values(SHELLS.iter()))]
    shell: String,
  },

  /// Estimates time trees from an initial tree topology, a set of date constraints (e.g. tip dates), and an alignment (optional).
  Timetree(TreetimeTimetreeArgs),

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
}

#[derive(Parser, Debug)]
pub struct TreetimeAncestralReassortmentGraphArgs;

pub fn generate_shell_completions(shell: &str) -> Result<(), Report> {
  let mut command = TreetimeArgs::command();

  if shell.to_lowercase() == "fig" {
    generate(Fig, &mut command, "treetime", &mut io::stdout());
    return Ok(());
  }

  let generator = <Shell as ArgEnum>::from_str(&shell.to_lowercase(), true)
    .map_err(|err| eyre!("{}: Possible values: {}", err, SHELLS.join(", ")))?;

  let bin_name = command.get_name().to_owned();

  generate(generator, &mut command, bin_name, &mut io::stdout());

  Ok(())
}

pub fn treetime_parse_cli_args() -> Result<TreetimeArgs, Report> {
  let args = TreetimeArgs::parse();

  // --verbosity=<level> and --silent take priority over -v and -q
  let filter_level = if args.silent {
    LevelFilter::Off
  } else {
    match args.verbosity {
      None => args.verbose.log_level_filter(),
      Some(verbosity) => verbosity,
    }
  };

  setup_logger(filter_level);

  Ok(args)
}
