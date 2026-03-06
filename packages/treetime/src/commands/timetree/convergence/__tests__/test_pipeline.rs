#[cfg(test)]
mod tests {
  use crate::commands::timetree::args::TreetimeTimetreeArgs;
  use crate::commands::timetree::run::run_timetree_estimation;
  use eyre::Report;
  use std::path::PathBuf;

  #[test]
  fn test_pipeline_timetree_convergence() -> Result<(), Report> {
    let root = project_root();
    let outdir = root.join("tmp/test-convergence-pipeline");
    std::fs::create_dir_all(&outdir)?;
    let tracelog_path = outdir.join("tracelog.csv");

    let args = TreetimeTimetreeArgs {
      input_fastas: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      dates: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      max_iter: 3,
      outdir: outdir.clone(),
      tracelog: Some(tracelog_path.clone()),
      ..TreetimeTimetreeArgs::default()
    };

    run_timetree_estimation(&args)?;

    // Verify tracelog was written and contains data
    let csv_content = std::fs::read_to_string(&tracelog_path)?;
    let lines: Vec<&str> = csv_content.lines().collect();
    assert!(lines.len() >= 2, "Tracelog must have header + at least 1 data row");

    // Verify header contains expected columns
    let header = lines[0];
    assert!(header.contains("n_diff"), "Tracelog header missing n_diff column");
    assert!(header.contains("n_resolved"), "Tracelog header missing n_resolved column");

    // Verify output tree files exist and are non-empty
    let nwk_path = outdir.join("timetree.nwk");
    let nex_path = outdir.join("timetree.nexus");
    assert!(nwk_path.exists(), "Output newick file must exist");
    assert!(nex_path.exists(), "Output nexus file must exist");
    assert!(
      std::fs::metadata(&nwk_path)?.len() > 0,
      "Output newick file must be non-empty"
    );

    Ok(())
  }

  fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .map(PathBuf::from)
      .expect("project has workspace root")
  }
}
