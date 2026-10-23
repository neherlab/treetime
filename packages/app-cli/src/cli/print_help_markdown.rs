use crate::cli::treetime_cli::TreetimeArgs;
use clap::{Command, CommandFactory};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use regex::{Captures, Regex};
use std::borrow::Cow;
use std::io::{self, Write};
use std::iter;

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
  let help = usage_code_blocks(&help)?;
  let help = examples_code_blocks(&help)?;

  Ok(help.into_owned())
}

fn usage_code_blocks(help: &str) -> Result<Cow<'_, str>, Report> {
  let usage = Regex::new(r"(?s)(## `([^`]+)`\n.*?)\*\*Usage:\*\* `([^`]*)`")?;
  Ok(usage.replace_all(help, |captures: &Captures<'_>| {
    let path = &captures[2];
    let lines = captures[3].lines().map(str::trim).collect_vec();
    let lines = match (path.rsplit_once(' '), lines.split_first()) {
      (Some((parent, _)), Some((&first, rest))) => {
        let first = first
          .strip_prefix(parent)
          .and_then(|line| line.strip_prefix(' '))
          .filter(|line| line.starts_with(path))
          .unwrap_or(first);
        iter::once(first).chain(rest.iter().copied()).join("\n")
      },
      _ => lines.join("\n"),
    };
    format!("{}**Usage:**\n\n```sh\n{lines}\n```", &captures[1])
  }))
}

fn examples_code_blocks(help: &str) -> Result<Cow<'_, str>, Report> {
  let examples = Regex::new(r"(?m)^Examples:\n((?:  .+\n?)+)")?;
  Ok(examples.replace_all(help, |captures: &Captures<'_>| {
    let lines = captures[1].lines().map(|line| line.trim_start_matches("  ")).join("\n");
    format!("**Examples:**\n\n```sh\n{lines}\n```\n")
  }))
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
  use crate::__tests__::test_cli_support::tests::helpers::{EVERY, find_command, treetime_command};
  use crate::cli::print_help_markdown::help_markdown;
  use crate::cli::treetime_cli::TreetimeArgs;
  use clap::CommandFactory;
  use eyre::Report;
  use helpers::{code_block, sections};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use regex::Regex;

  #[test]
  fn test_print_help_markdown_usage_names_the_command_path_once() -> Result<(), Report> {
    let help = help_markdown()?;
    let expected_misnamed: Vec<(String, String)> = vec![];
    let expected_sampled = vec![
      (
        "treetime timetree",
        "treetime timetree --tree <FILE> --metadata <FILE> [(--alignment <FILE>)... | --sequence-length <N>] (-O <DIR> | --output-* <FILE>) [OPTIONS]",
      ),
      ("treetime timetree", "treetime timetree --config <FILE> [OPTIONS]"),
      ("treetime schema write", "treetime schema write [--for <TARGET>] [-o <PATH>] [OPTIONS]"),
      ("treetime help-markdown", "treetime help-markdown [OPTIONS] > <FILE>"),
    ];

    let sections = sections(&help);
    let usage = sections
      .iter()
      .flat_map(|(path, section)| {
        code_block(section, "**Usage:**")
          .unwrap_or_default()
          .lines()
          .map(|line| (path.clone(), line.to_owned()))
          .collect_vec()
      })
      .collect_vec();
    let misnamed = usage
      .iter()
      .filter(|(path, line)| line != path && !line.starts_with(&format!("{path} ")))
      .cloned()
      .collect_vec();
    let without_usage = sections
      .iter()
      .filter(|(_, section)| code_block(section, "**Usage:**").is_none())
      .count();
    let sampled = usage
      .iter()
      .filter(|(path, _)| {
        ["treetime timetree", "treetime schema write", "treetime help-markdown"].contains(&path.as_str())
      })
      .map(|(path, line)| (path.as_str(), line.as_str()))
      .collect_vec();

    assert_eq!(
      (expected_misnamed, 0, expected_sampled),
      (misnamed, without_usage, sampled)
    );
    Ok(())
  }

  #[test]
  fn test_print_help_markdown_examples_are_code_blocks() -> Result<(), Report> {
    let root = treetime_command();
    let help = help_markdown()?;
    let expected = EVERY
      .iter()
      .map(|path| {
        let after_long_help = find_command(&root, path)
          .unwrap()
          .get_after_long_help()
          .map(ToString::to_string)
          .unwrap_or_default();
        let lines = after_long_help
          .lines()
          .filter_map(|line| line.strip_prefix("  "))
          .join("\n");
        (format!("treetime {path}"), lines)
      })
      .sorted()
      .collect_vec();

    let actual = sections(&help)
      .into_iter()
      .filter_map(|(path, section)| code_block(&section, "**Examples:**").map(|block| (path, block)))
      .sorted()
      .collect_vec();

    assert_eq!((expected, 0), (actual, help.matches("\nExamples:\n").count()));
    Ok(())
  }

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
        let jobs_help = help.split("--jobs <N>").nth(1).unwrap_or_default();
        let jobs_help = next_argument.split(jobs_help).next().unwrap_or_default();
        (command.get_name().to_owned(), jobs_help.contains("default:"))
      })
      .collect_vec();

    let help = help_markdown()?;
    let markdown_defaults = help
      .split("`--jobs <N>`")
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
          ("ancestral", false),
          ("clock", false),
          ("mugration", false),
          ("homoplasy", false),
          ("optimize", false),
          ("prune", false),
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

  mod helpers {
    use itertools::Itertools;

    pub(super) fn sections(help: &str) -> Vec<(String, String)> {
      help
        .split("\n## `")
        .skip(1)
        .filter_map(|section| {
          section
            .split_once('`')
            .map(|(path, rest)| (path.to_owned(), rest.to_owned()))
        })
        .collect_vec()
    }

    pub(super) fn code_block(section: &str, label: &str) -> Option<String> {
      let (_, block) = section.split_once(&format!("{label}\n\n```sh\n"))?;
      block.split_once("\n```").map(|(block, _)| block.to_owned())
    }
  }
}
