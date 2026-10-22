use crate::cli::color::{ColorArgs, color_choice_of_args};
use crate::cli::pipeline::runner::run_pipeline_command;
use crate::cli::print_debug_info::print_debug_info;
use crate::cli::print_help_markdown::print_help_markdown;
use crate::cli::progress::{BarProgress, TextProgress};
use crate::cli::rtt_chart::print_clock_regression_chart;
use crate::cli::schema::generate_schema;
use crate::cli::treetime_cli::{
  CliParseError, ParsedCli, TreetimeCommands, TreetimeSchemaArgs, TreetimeSchemaCommands, generate_shell_completions,
  resolve_args, treetime_parse_cli_args,
};
use crate::cli::verbosity::Verbosity;
use app_commands::command::{CommandArgs, RawConfig};
use app_commands::commands::clock::run::run_clock;
use app_commands::commands::shared::argument_error::ArgumentError;
use app_commands::config::source::InvalidConfig;
use eyre::Report;
use itertools::Itertools;
use log::{debug, info};
use std::env;
use std::process::ExitCode;
use std::sync::Arc;
use treetime::cancel::NoopCancel;
use treetime::progress::{LogSink, NoopProgress, StageSink};
use treetime_utils::error::report_to_string;
use treetime_utils::init::clap_styles::styles;
use treetime_utils::init::global::{install_error_hook, install_panic_hook, setup_logger};
use treetime_utils::init::thread_pool::init_thread_pool;
use treetime_utils::io::console::is_tty;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::o;

const USAGE_ERROR_EXIT_CODE: u8 = 2;

pub fn run_cli() -> ExitCode {
  let argv = env::args_os().collect_vec();
  let color = color_choice_of_args(&argv);
  let parsed = match treetime_parse_cli_args(&argv, color) {
    Ok(parsed) => parsed,
    Err(CliParseError::Clap(err)) => err.exit(),
    Err(CliParseError::Argument(report)) => return report_failure(&report, ColorArgs { color }.stderr_has_colors()),
  };
  let colors = parsed.args.color.stderr_has_colors();
  install_panic_hook(colors);
  match install_error_hook(colors).and_then(|()| run_parsed(parsed, colors)) {
    Ok(()) => ExitCode::SUCCESS,
    Err(report) => report_failure(&report, colors),
  }
}

fn run_parsed(parsed: ParsedCli, colors: bool) -> Result<(), Report> {
  setup_logger(parsed.args.verbosity.get_filter_level(), colors);
  for message in &parsed.v0_messages {
    debug!("{message}");
  }
  let args = resolve_args(parsed)?;

  info!("# Command line arguments");
  info!("{}", json_write_str(&args, JsonPretty(true))?);

  if let Some(jobs) = args.command.jobs() {
    init_thread_pool(jobs)?;
  }

  let (stages, log) = make_progress(&args.verbosity, colors)?;
  run_command(args.command, &*stages, &*log)
}

#[cfg_attr(
  dylint_lib = "custom",
  expect(debug_remnants, reason = "the CLI reports a failed run on standard error")
)]
fn report_failure(report: &Report, colors: bool) -> ExitCode {
  if let Some(message) = usage_error(report) {
    let style = *styles().get_error();
    let label = if colors {
      format!("{}error:{}", style.render(), style.render_reset())
    } else {
      o!("error:")
    };
    eprintln!("{label} {message}");
    ExitCode::from(USAGE_ERROR_EXIT_CODE)
  } else {
    eprintln!("Error: {report:?}");
    ExitCode::FAILURE
  }
}

fn usage_error(report: &Report) -> Option<String> {
  let message = report_to_string(report);
  if let Some(invalid) = report.chain().find_map(|error| error.downcast_ref::<InvalidConfig>()) {
    return Some(format!("{message}\n\n{}", invalid.rendered));
  }
  report
    .chain()
    .any(|error| error.is::<ArgumentError>())
    .then_some(message)
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
      let clock_args = clock_args.into_args()?;
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
    TreetimeCommands::Schema {
      command: TreetimeSchemaCommands::Write(TreetimeSchemaArgs { target, output }),
    } => {
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

fn make_progress(verbosity: &Verbosity, colors: bool) -> Result<(Arc<dyn StageSink>, Arc<dyn LogSink>), Report> {
  Ok(match verbosity.get_log_level() {
    None => shared_sink(NoopProgress),
    Some(min_level) => {
      if !verbosity.no_progress && is_tty() {
        shared_sink(BarProgress::new(min_level, colors)?)
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
