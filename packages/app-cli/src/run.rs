use crate::cli::examples::run_examples_get;
use crate::cli::pipeline::runner::run_pipeline_command;
use crate::cli::print_debug_info::print_debug_info;
use crate::cli::print_help_markdown::print_help_markdown;
use crate::cli::progress::{BarProgress, TextProgress};
use crate::cli::rtt_chart::print_clock_regression_chart;
use crate::cli::schema::generate_schema;
use crate::cli::treetime_cli::{
  TreetimeCommands, TreetimeExamplesCommands, TreetimeSchemaArgs, generate_shell_completions, treetime_parse_cli_args,
};
use crate::cli::verbosity::Verbosity;
use app_commands::command::CommandArgs;
use app_commands::commands::clock::args::TreetimeClockArgs;
use app_commands::commands::clock::run::run_clock;
use eyre::Report;
use log::info;
use std::env;
use std::sync::Arc;
use treetime::cancel::NoopCancel;
use treetime::progress::{LogSink, NoopProgress, StageSink};
use treetime_utils::init::global::setup_logger;
use treetime_utils::init::thread_pool::init_thread_pool;
use treetime_utils::io::console::is_tty;
use treetime_utils::io::json::{JsonPretty, json_write_str};

pub fn run_cli() -> Result<(), Report> {
  let args = treetime_parse_cli_args(env::args_os())?;
  setup_logger(args.verbosity.get_filter_level());

  info!("# Command line arguments");
  info!("{}", json_write_str(&args, JsonPretty(true))?);

  init_thread_pool(args.jobs.jobs)?;

  let (stages, log) = make_progress(&args.verbosity)?;
  run_command(args.command, &*stages, &*log)
}

pub(crate) fn run_command(command: TreetimeCommands, stages: &dyn StageSink, log: &dyn LogSink) -> Result<(), Report> {
  match command {
    TreetimeCommands::Completions { shell } => {
      generate_shell_completions(&shell)?;
    },
    TreetimeCommands::Timetree(args) => CommandArgs::try_from(*args)?.execute(&NoopCancel, stages, log)?,
    TreetimeCommands::Optimize(args) => CommandArgs::try_from(args)?.execute(&NoopCancel, stages, log)?,
    TreetimeCommands::Prune(args) => CommandArgs::try_from(args)?.execute(&NoopCancel, stages, log)?,
    TreetimeCommands::Ancestral(args) => CommandArgs::try_from(args)?.execute(&NoopCancel, stages, log)?,
    TreetimeCommands::Clock(clock_args) => {
      let clock_args = TreetimeClockArgs::try_from(clock_args)?;
      let result = run_clock(&clock_args, &NoopCancel, stages, log)?;
      if is_tty() {
        print_clock_regression_chart(&result.regression_results, &result.clock_model)?;
      }
    },
    TreetimeCommands::Homoplasy(args) => CommandArgs::try_from(args)?.execute(&NoopCancel, stages, log)?,
    TreetimeCommands::Mugration(args) => CommandArgs::try_from(args)?.execute(&NoopCancel, stages, log)?,
    TreetimeCommands::Pipeline(pipeline_args) => {
      run_pipeline_command(&pipeline_args, stages, log)?;
    },
    TreetimeCommands::Examples {
      command: TreetimeExamplesCommands::Get(examples_args),
    } => {
      run_examples_get(&examples_args)?;
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

fn make_progress(verbosity: &Verbosity) -> Result<(Arc<dyn StageSink>, Arc<dyn LogSink>), Report> {
  Ok(match verbosity.get_log_level() {
    None => shared_sink(NoopProgress),
    Some(min_level) => {
      if !verbosity.no_progress && is_tty() {
        shared_sink(BarProgress::new(min_level)?)
      } else {
        shared_sink(TextProgress::new(min_level))
      }
    },
  })
}

fn shared_sink<S: StageSink + LogSink + 'static>(sink: S) -> (Arc<dyn StageSink>, Arc<dyn LogSink>) {
  let sink = Arc::new(sink);
  let stages: Arc<dyn StageSink> = Arc::<S>::clone(&sink);
  (stages, sink)
}
