use crate::io::fs::filename_maybe;
use crate::utils::datetime::{date_format_precise, date_now};
use color_eyre::owo_colors::{OwoColorize, Style};
use env_logger::Env;
use log::{Level, LevelFilter, Record};
use std::env;
use std::io::Write;

fn get_current_exe_filename() -> Option<String> {
  env::current_exe().ok().and_then(filename_maybe)
}

fn get_file_line(record: &Record) -> String {
  let file = record.file().and_then(filename_maybe);
  match (file, record.line()) {
    (Some(file), None) => format!("{file}:",),
    (Some(file), Some(line)) => format!("{file}:{line:}:"),
    _ => "".to_owned(),
  }
  .dimmed()
  .to_string()
}

fn log_level_str(record: &Record) -> String {
  let mut level_str = record.level().to_string();
  level_str.truncate(1);
  level_str
}

fn color_log_level(record: &Record) -> String {
  let level_str = match record.level() {
    Level::Error => log_level_str(record).red().to_string(),
    Level::Warn => log_level_str(record).yellow().to_string(),
    Level::Info => log_level_str(record).cyan().dimmed().to_string(),
    Level::Debug => log_level_str(record).green().dimmed().to_string(),
    Level::Trace => log_level_str(record).dimmed().to_string(),
  };
  format!("{:}{level_str}{:}", "[".dimmed(), "]".dimmed())
}

pub fn setup_logger(filter_level: LevelFilter) {
  env_logger::Builder::from_env(Env::default().default_filter_or("warn"))
    .filter_level(filter_level)
    .format(|buf, record| {
      let current_exe = get_current_exe_filename().unwrap_or_default().dimmed().to_string();
      let file_line = get_file_line(record);
      let level = color_log_level(record);
      let date = date_format_precise(&date_now()).dimmed().to_string();
      let args = record.args();
      writeln!(buf, "{date} {level:} {file_line:} {args}")?;
      Ok(())
    })
    .init();
}

pub fn global_init() {
  color_eyre::config::HookBuilder::default()
    .theme(
      color_eyre::config::Theme::dark()
        .dependency_code(Style::new().dimmed())
        .file(Style::new().green())
        .line_number(Style::new().yellow())
        .panic_file(Style::new().green())
        .panic_line_number(Style::new().yellow())
        .panic_message(Style::new().bright_red().bold())
        .active_line(Style::new().cyan())
        .hidden_frames(Style::new().dimmed())
        .code_hash(Style::new().hidden()),
    )
    .panic_section(format!(
      "If you think it's a bug, consider reporting at: '{}/issues'",
      env!("CARGO_PKG_REPOSITORY"),
    ))
    .add_frame_filter(Box::new(|frames| {
      frames.retain(|frame| {
        let should_show_name = frame.name.as_ref().map_or(false, |name| {
          !HIDDEN_CRATE_NAME_PREFIXES
            .iter()
            .any(|&prefix| name.starts_with(prefix) || name.starts_with(&format!("<{prefix}")))
        });

        let should_show_file = !frame.filename.as_ref().map_or(false, |filename| {
          HIDDEN_CRATE_PATH_PREFIXES
            .iter()
            .any(|&prefix| filename.starts_with(prefix))
        });

        should_show_file && should_show_name
      });
    }))
    .install()
    .expect("color_eyre initialization failed");
}

const HIDDEN_CRATE_NAME_PREFIXES: &[&str] = &[
  "__rust_try",
  "alloc::",
  "color_eyre::",
  "core::",
  "crossbeam::",
  "eyre::",
  "rayon::",
  "rayon_core::",
  "rustc::",
  "std::",
  "tokio::",
];

const HIDDEN_CRATE_PATH_PREFIXES: &[&str] = &["/rustc/"];
