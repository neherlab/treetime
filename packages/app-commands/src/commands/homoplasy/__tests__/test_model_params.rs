#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{
    HKY85_PARAMS, RecordingLog, WrittenGtr, assert_written_gtr, hky85_model, project_root,
  };
  use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
  use crate::commands::homoplasy::run::run_homoplasy;
  use clap::Parser;
  use eyre::Report;
  use itertools::Itertools;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::{LogLevel, NoopProgress};
  use treetime_utils::io::json::json_read_str;

  #[test]
  fn test_homoplasy_model_params_reach_the_model() -> Result<(), Report> {
    let out = tempdir()?;
    let data = project_root().join("data/flu/h3n2/20");
    let tree = format!("--tree={}", data.join("tree.nwk").display());
    let alignment = format!("--alignment={}", data.join("aln.fasta.xz").display());
    let output = format!("--output-homoplasy-stats={}", out.path().join("stats.json").display());
    let mut argv = vec![
      "homoplasy",
      &tree,
      &alignment,
      "--method-anc=marginal",
      "--model=hky85",
      "--model-params",
    ];
    argv.extend(HKY85_PARAMS);
    argv.push(&output);
    let args = TreetimeHomoplasyArgs::try_from(TreetimeHomoplasyArgsRaw::try_parse_from(argv)?)?;
    let log = RecordingLog::default();

    run_homoplasy(&args, &NoopCancel, &NoopProgress, &log)?;

    let models: Vec<WrittenGtr> = log
      .messages(LogLevel::Info)
      .iter()
      .filter_map(|message| message.strip_prefix("GTR model initialized:\n"))
      .map(json_read_str)
      .try_collect()?;
    let [model] = models.as_slice() else {
      panic!("one alignment gives one model");
    };
    assert_written_gtr(&hky85_model(), model)
  }
}
