use crate::cli::treetime_cli::TreetimeArgs;
use clap::{Command, CommandFactory};
use eyre::{Report, WrapErr};
use regex::Regex;
use std::borrow::Cow;
use std::io::{self, Write};

pub(crate) fn print_help_markdown() -> Result<(), Report> {
  let help = help_markdown()?;

  writeln!(io::stdout().lock(), "{help}").wrap_err("When writing help to standard output")?;
  Ok(())
}

fn help_markdown() -> Result<String, Report> {
  let command = markdown_command();
  let help = clap_markdown::help_markdown_command(&command);

  let help = replace(&help, "# Command-Line Help for `treetime`", "")?;

  let help = replace(
    &help,
    "This document contains the help content for the `treetime` command-line program.",
    r#"
This document contains the automatically generated reference documentation for command-line arguments of the latest version of TreeTime CLI.

If you have TreeTime CLI installed, you can type `treetime --help` to read the latest documentation for your installed version of TreeTime. To generate this document in markdown format, run `treetime help-markdown > reference.md`
  "#,
  )?;

  let help = replace(&help, r"(.*)\x{2014} REMOVED(.*)", "")?;
  let help = replace(&help, r"(.*)\x{2014} RENAMED(.*)", "")?;

  Ok(help.into_owned())
}

fn markdown_command() -> Command {
  hide_jobs_default_recursive(TreetimeArgs::command())
}

fn hide_jobs_default_recursive(command: Command) -> Command {
  command
    .mut_args(|arg| {
      if arg.get_id() == "jobs" {
        arg.default_value(None::<&str>)
      } else {
        arg
      }
    })
    .mut_subcommands(hide_jobs_default_recursive)
}

fn replace<'t>(text: &'t str, what: &str, with_what: &str) -> Result<Cow<'t, str>, Report> {
  let res = Regex::new(what)
    .wrap_err_with(|| format!("When compiling regex: {what}"))?
    .replace_all(text, with_what);
  Ok(res)
}

#[cfg(test)]
mod tests {
  use crate::cli::print_help_markdown::help_markdown;
  use crate::cli::treetime_cli::TreetimeArgs;
  use clap::CommandFactory;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use regex::Regex;

  #[test]
  fn test_print_help_markdown_jobs_default_is_hidden() -> Result<(), Report> {
    let next_argument = Regex::new(r"\n +-")?;
    let mut root = TreetimeArgs::command();
    root.build();
    let terminal_defaults = root
      .get_subcommands_mut()
      .filter(|command| command.get_arguments().any(|arg| arg.get_id() == "jobs"))
      .map(|command| {
        let help = command.render_long_help().to_string();
        let jobs_help = help.split("--jobs <JOBS>").nth(1).unwrap_or_default();
        let jobs_help = next_argument.split(jobs_help).next().unwrap_or_default();
        (command.get_name().to_owned(), jobs_help.contains("default:"))
      })
      .collect_vec();

    let help = help_markdown()?;
    let markdown_defaults = help
      .split("`--jobs <JOBS>`")
      .skip(1)
      .map(|jobs_help| {
        let jobs_help = jobs_help.split("\n* `").next().unwrap_or_default();
        jobs_help
          .split("\n## ")
          .next()
          .unwrap_or_default()
          .contains("Default value:")
      })
      .collect_vec();

    assert_eq!(
      (
        vec![
          ("timetree", false),
          ("optimize", false),
          ("prune", false),
          ("ancestral", false),
          ("clock", false),
          ("homoplasy", false),
          ("mugration", false),
          ("pipeline", false),
        ],
        vec![false; 8]
      ),
      (
        terminal_defaults
          .iter()
          .map(|(name, shown)| (name.as_str(), *shown))
          .collect_vec(),
        markdown_defaults
      )
    );
    Ok(())
  }
}
