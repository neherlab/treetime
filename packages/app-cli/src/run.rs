use crate::cli::pipeline::runner::run_pipeline_command;
use crate::cli::print_debug_info::print_debug_info;
use crate::cli::print_help_markdown::print_help_markdown;
use crate::cli::progress::{BarProgress, TextProgress};
use crate::cli::rtt_chart::print_clock_regression_chart;
use crate::cli::schema::generate_schema;
use crate::cli::treetime_cli::{
  TreetimeCommands, TreetimeSchemaArgs, generate_shell_completions, treetime_parse_cli_args,
};
use crate::cli::verbosity::Verbosity;
use crate::commands::ancestral::args::TreetimeAncestralArgs;
use crate::commands::ancestral::run::run_ancestral_reconstruction;
use crate::commands::clock::args::TreetimeClockArgs;
use crate::commands::clock::run::run_clock;
use crate::commands::homoplasy::args::TreetimeHomoplasyArgs;
use crate::commands::homoplasy::run::run_homoplasy;
use crate::commands::mugration::args::TreetimeMugrationArgs;
use crate::commands::mugration::run::run_mugration;
use crate::commands::optimize::args::TreetimeOptimizeArgs;
use crate::commands::optimize::run::run_optimize;
use crate::commands::prune::args::TreetimePruneArgs;
use crate::commands::prune::run::run_prune;
use crate::commands::timetree::args::TreetimeTimetreeArgs;
use crate::commands::timetree::run::run_timetree_estimation;
use eyre::Report;
use log::info;
use treetime::cancel::NoopCancel;
use treetime::progress::{NoopProgress, ProgressSink};
use treetime_utils::io::console::is_tty;
use treetime_utils::io::json::{JsonPretty, json_write_str};

pub fn run_cli() -> Result<(), Report> {
  let args = treetime_parse_cli_args()?;

  info!("# Command line arguments");
  info!("{}", json_write_str(&args, JsonPretty(true))?);

  if args.jobs.jobs == 1 {
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .use_current_thread()
      .build_global()?;
  } else {
    rayon::ThreadPoolBuilder::new()
      .num_threads(args.jobs.jobs)
      .build_global()?;
  }

  let progress = make_progress(&args.verbosity)?;

  match args.command {
    TreetimeCommands::Completions { shell } => {
      generate_shell_completions(&shell)?;
    },
    TreetimeCommands::Timetree(timetree_args) => {
      let timetree_args = TreetimeTimetreeArgs::try_from(*timetree_args)?;
      run_timetree_estimation(&timetree_args, &NoopCancel, &*progress)?;
    },
    TreetimeCommands::Optimize(optimize_args) => {
      let optimize_args = TreetimeOptimizeArgs::try_from(optimize_args)?;
      run_optimize(&optimize_args, &NoopCancel, &*progress)?;
    },
    TreetimeCommands::Prune(prune_args) => {
      let prune_args = TreetimePruneArgs::try_from(prune_args)?;
      run_prune(&prune_args, &NoopCancel, &*progress)?;
    },
    TreetimeCommands::Ancestral(ancestral_args) => {
      let ancestral_args = TreetimeAncestralArgs::try_from(ancestral_args)?;
      run_ancestral_reconstruction(&ancestral_args, &NoopCancel, &*progress)?;
    },
    TreetimeCommands::Clock(clock_args) => {
      let clock_args = TreetimeClockArgs::try_from(clock_args)?;
      let result = run_clock(&clock_args, &NoopCancel, &*progress)?;
      if is_tty() {
        print_clock_regression_chart(&result.regression_results, &result.clock_model)?;
      }
    },
    TreetimeCommands::Homoplasy(homoplasy_args) => {
      let homoplasy_args = TreetimeHomoplasyArgs::try_from(homoplasy_args)?;
      run_homoplasy(&homoplasy_args, &NoopCancel, &*progress)?;
    },
    TreetimeCommands::Mugration(mugration_args) => {
      let mugration_args = TreetimeMugrationArgs::try_from(mugration_args)?;
      run_mugration(&mugration_args, &NoopCancel, &*progress)?;
    },
    TreetimeCommands::Pipeline(pipeline_args) => {
      run_pipeline_command(&pipeline_args, &*progress)?;
    },
    TreetimeCommands::Arg(_) => {},
    TreetimeCommands::Schema(TreetimeSchemaArgs { target, output }) => {
      generate_schema(target, output.as_ref())?;
    },
    TreetimeCommands::HelpMarkdown => {
      print_help_markdown()?;
    },
    TreetimeCommands::Debug => {
      print_debug_info()?;
    },
  }

  Ok(())
}

fn make_progress(verbosity: &Verbosity) -> Result<Box<dyn ProgressSink>, Report> {
  Ok(match verbosity.get_log_level() {
    None => Box::new(NoopProgress),
    Some(min_level) => {
      if !verbosity.no_progress && is_tty() {
        Box::new(BarProgress::new(min_level)?)
      } else {
        Box::new(TextProgress::new(min_level))
      }
    },
  })
}
