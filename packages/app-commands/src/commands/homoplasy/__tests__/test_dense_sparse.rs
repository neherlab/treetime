#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::run;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rstest]
  #[case::zika_86("zika/86")]
  #[case::flu_h3n2_20("flu/h3n2/20")]
  #[trace]
  fn test_dense_sparse_homoplasy_statistics_agree(#[case] dataset: &str) -> Result<(), Report> {
    assert_eq!(run(dataset, true)?, run(dataset, false)?);
    Ok(())
  }

  #[ignore = "dense and sparse reconstruct different internal residues on rsv/a/20, see kb/issues/M-ancestral-sparse-dense-internal-residues-diverge.md"]
  #[test]
  fn test_dense_sparse_homoplasy_statistics_agree_with_gaps_and_unknown_characters() -> Result<(), Report> {
    assert_eq!(run("rsv/a/20", true)?, run("rsv/a/20", false)?);
    Ok(())
  }

  mod helpers {
    use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
    use crate::commands::homoplasy::result::HomoplasyResult;
    use crate::commands::homoplasy::run::run_homoplasy;
    use clap::Parser;
    use eyre::Report;
    use std::path::Path;
    use tempfile::tempdir;
    use treetime::cancel::NoopCancel;
    use treetime::progress::NoopProgress;

    pub(super) fn run(dataset: &str, dense: bool) -> Result<HomoplasyResult, Report> {
      let out = tempdir()?;
      let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(dataset);
      let raw = TreetimeHomoplasyArgsRaw::try_parse_from([
        "homoplasy",
        &format!("--tree={}", data.join("tree.nwk").display()),
        &format!("--alignment={}", data.join("aln.fasta.xz").display()),
        "--model=jc69",
        "--method-anc=marginal",
        &format!("--dense={dense}"),
        &format!("--output-homoplasy-stats={}", out.path().join("stats.json").display()),
      ])?;
      let args = TreetimeHomoplasyArgs::try_from(raw)?;
      run_homoplasy(&args, &NoopCancel, &NoopProgress, &NoopProgress)
    }
  }
}
