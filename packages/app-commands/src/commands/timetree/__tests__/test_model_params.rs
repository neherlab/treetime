#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{
    HKY85_PARAMS, WrittenModels, assert_written_gtr, hky85_model, project_root,
  };
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::run::run_timetree_estimation;
  use clap::Parser;
  use eyre::Report;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::io::json::json_read_file;

  #[test]
  fn test_timetree_model_params_reach_the_model() -> Result<(), Report> {
    let out = tempdir()?;
    let data = project_root().join("data/flu/h3n2/20");
    let model_file = out.path().join("model.json");
    let tree = format!("--tree={}", data.join("tree.nwk").display());
    let alignment = format!("--alignment={}", data.join("aln.fasta.xz").display());
    let metadata = format!("--metadata={}", data.join("metadata.tsv").display());
    let output = format!("--output-model={}", model_file.display());
    let mut argv = vec![
      "timetree",
      &tree,
      &alignment,
      &metadata,
      "--branch-length-mode=marginal",
      "--max-iter=1",
      "--seed=1",
      "--model=hky85",
      "--model-params",
    ];
    argv.extend(HKY85_PARAMS);
    argv.push(&output);
    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw::try_parse_from(argv)?)?;

    run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let written: WrittenModels = json_read_file(&model_file)?;
    let [model] = written.models.as_slice() else {
      panic!("one alignment gives one model");
    };
    assert_written_gtr(&hky85_model(), model)
  }
}
