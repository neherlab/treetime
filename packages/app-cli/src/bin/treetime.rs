<<<<<<< HEAD
<<<<<<< HEAD
#[cfg(any(
  all(target_arch = "x86_64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "x86_64", target_os = "linux", target_env = "musl"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "gnu"),
  all(target_arch = "aarch64", target_os = "linux", target_env = "musl"),
))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

=======
use app_cli::cli::progress::BarProgress;
>>>>>>> 04a8f9b2 (feat(app-cli): add indicatif progress bar for TTY, NoopProgress for pipes)
use app_cli::cli::rtt_chart::{
  print_clock_regression_chart, write_clock_regression_chart_png, write_clock_regression_chart_svg,
};
<<<<<<< HEAD
use app_cli::cli::treetime_cli::{TreetimeCommands, TreetimeSchemaArgs, generate_shell_completions, treetime_parse_cli_args};
=======
use app_cli::cli::progress::{BarProgress, TextProgress};
use app_cli::cli::rtt_chart::{
  print_clock_regression_chart, write_clock_regression_chart_png, write_clock_regression_chart_svg,
};
use app_cli::cli::treetime_cli::{
  TreetimeCommands, TreetimeSchemaArgs, generate_shell_completions, treetime_parse_cli_args,
};
use app_cli::cli::verbosity::Verbosity;
>>>>>>> d9cec2f3 (feat(cli): add --no-progress flag, BarProgress with suspend, TextProgress)
=======
use app_cli::cli::treetime_cli::{
  TreetimeCommands, TreetimeSchemaArgs, generate_shell_completions, treetime_parse_cli_args,
};
>>>>>>> f1127239 (refactor: format)
use ctor::ctor;
use eyre::Report;
use log::info;
use treetime::commands::ancestral::run::run_ancestral_reconstruction;
use treetime::commands::clock::run::run_clock;
use treetime::commands::homoplasy::run::run_homoplasy;
use treetime::commands::mugration::run::run_mugration;
use treetime::commands::optimize::run::run_optimize;
use treetime::commands::prune::run::run_prune;
use treetime::commands::timetree::run::run_timetree_estimation;
<<<<<<< HEAD
use treetime::progress::{NoopProgress, ProgressSink};
=======
>>>>>>> f1127239 (refactor: format)
use treetime::schema::generate_schema;
use treetime_utils::init::global::global_init;
use treetime_utils::init::openblas::get_openblas_info_str;
use treetime_utils::io::console::is_tty;
use treetime_utils::io::json::{JsonPretty, json_write_str};

#[ctor]
fn init() {
  global_init();
}

fn make_progress(verbosity: &Verbosity) -> Box<dyn ProgressSink> {
  match verbosity.get_log_level() {
    None => Box::new(NoopProgress),
    Some(min_level) => {
      if !verbosity.no_progress && is_tty() {
        Box::new(BarProgress::new(min_level))
      } else {
        Box::new(TextProgress::new(min_level))
      }
    },
  }
}

fn main() -> Result<(), Report> {
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

  let progress = make_progress(&args.verbosity);

  match args.command {
    TreetimeCommands::Timetree(timetree_args) => {
      run_timetree_estimation(&timetree_args, &*progress)?;
    },
    TreetimeCommands::Optimize(optimize_args) => {
      run_optimize(&optimize_args, &*progress)?;
    },
    TreetimeCommands::Prune(prune_args) => {
      run_prune(&prune_args, &*progress)?;
    },
    TreetimeCommands::Ancestral(ancestral_args) => {
      run_ancestral_reconstruction(&ancestral_args, &*progress)?;
    },
    TreetimeCommands::Clock(clock_args) => {
      let result = run_clock(&clock_args, &*progress)?;
      let outdir = &clock_args.outdir;
      write_clock_regression_chart_svg(
        &result.regression_results,
        &result.clock_model,
        outdir.join("clock.svg"),
      )?;
      write_clock_regression_chart_png(
        &result.regression_results,
        &result.clock_model,
        outdir.join("clock.png"),
      )?;
      if is_tty() {
        print_clock_regression_chart(&result.regression_results, &result.clock_model)?;
      }
    },
    TreetimeCommands::Homoplasy(homoplasy_args) => {
      run_homoplasy(homoplasy_args)?;
    },
    TreetimeCommands::Mugration(mugration_args) => {
      run_mugration(&mugration_args, &*progress)?;
    },
    TreetimeCommands::Completions { shell } => {
      generate_shell_completions(&shell)?;
    },
    TreetimeCommands::Schema(TreetimeSchemaArgs { for_format, output }) => {
      generate_schema(&for_format, output.as_ref())?;
    },
    TreetimeCommands::Arg(arg_args) => {},
    TreetimeCommands::Debug => {
      println!("{}", get_openblas_info_str());
    },
  }

  Ok(())
}
