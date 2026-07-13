use crate::cli::treetime_cli::TreetimeArgs;
use clap::{Command, CommandFactory};
use eyre::{Report, WrapErr};
use regex::Regex;
use std::borrow::Cow;

pub fn print_help_markdown() -> Result<(), Report> {
  let help = help_markdown()?;

  println!("{help}");
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
  // clap-markdown 0.1.5 ignores Clap's hide_default_value setting.
  TreetimeArgs::command().mut_args(|arg| {
    if arg.get_id() == "jobs" {
      arg.default_value(None::<&str>)
    } else {
      arg
    }
  })
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

  #[test]
  fn test_print_help_markdown_jobs_default_is_hidden() -> Result<(), Report> {
    let mut command = TreetimeArgs::command();
    let help = command.render_long_help().to_string();
    let jobs_help = help
      .split("--jobs")
      .nth(1)
      .expect("jobs argument should appear in terminal help")
      .split("--verbosity")
      .next()
      .expect("verbosity argument should follow jobs in terminal help");
    assert!(!jobs_help.contains("default:"));

    let help = help_markdown()?;
    let jobs_help = help
      .split("`--jobs <JOBS>`")
      .nth(1)
      .expect("jobs argument should appear in Markdown help")
      .split("`--verbosity")
      .next()
      .expect("verbosity argument should follow jobs in Markdown help");
    assert!(!jobs_help.contains("Default value:"));
    Ok(())
  }
}
