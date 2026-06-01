use crate::cli::treetime_cli::TreetimeArgs;
use eyre::{Report, WrapErr};
use regex::Regex;
use std::borrow::Cow;

pub fn print_help_markdown() -> Result<(), Report> {
  let help = clap_markdown::help_markdown::<TreetimeArgs>();

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

  println!("{help}");
  Ok(())
}

fn replace<'t>(text: &'t str, what: &str, with_what: &str) -> Result<Cow<'t, str>, Report> {
  let res = Regex::new(what)
    .wrap_err_with(|| format!("When compiling regex: {what}"))?
    .replace_all(text, with_what);
  Ok(res)
}
