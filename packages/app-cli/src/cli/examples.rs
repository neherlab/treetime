use crate::cli::treetime_cli::TreetimeExamplesGetArgs;
use app_commands::examples_download::{DownloadProgress, download_examples, examples_url};
use eyre::{Report, WrapErr};
use indicatif::{ProgressBar, ProgressStyle};
use treetime_utils::io::console::is_tty;

const BAR_TEMPLATE: &str = "{spinner:.green} [{bar:30.cyan/dim}] {bytes}/{total_bytes} {msg}";

pub(crate) fn run_examples_get(args: &TreetimeExamplesGetArgs) -> Result<(), Report> {
  let url = args.url.clone().unwrap_or_else(examples_url);
  let bar = is_tty().then(download_bar).transpose()?;
  let result = download_examples(&url, &args.output_dir, |DownloadProgress { received, total }| {
    if let Some(bar) = &bar {
      if let Some(total) = total {
        bar.set_length(total);
      }
      bar.set_position(received);
    }
  });
  if let Some(bar) = &bar {
    bar.finish_and_clear();
  }
  result
}

fn download_bar() -> Result<ProgressBar, Report> {
  let bar = ProgressBar::new(0);
  bar.set_style(
    ProgressStyle::with_template(BAR_TEMPLATE)
      .wrap_err("When parsing the progress bar template")?
      .progress_chars("=> "),
  );
  bar.set_message("example datasets");
  Ok(bar)
}
