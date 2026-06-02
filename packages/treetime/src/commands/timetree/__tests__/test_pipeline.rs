#[cfg(test)]
mod tests {
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output::OutputArgs;
  use crate::commands::timetree::args::TreetimeTimetreeArgs;
  use crate::commands::timetree::run::run_timetree_estimation;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use std::fs::read_to_string;
  use std::path::PathBuf;

  #[test]
  fn test_pipeline_timetree_convergence() -> Result<(), Report> {
    let root = project_root();
    let outdir = root.join("tmp/test-convergence-pipeline");
    std::fs::create_dir_all(&outdir)?;
    let tracelog_path = outdir.join("tracelog.csv");

    let args = TreetimeTimetreeArgs {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      max_iter: 3,
      output: OutputArgs {
        outdir: outdir.clone(),
        ..Default::default()
      },
      tracelog: Some(tracelog_path.clone()),
      ..TreetimeTimetreeArgs::default()
    };

    run_timetree_estimation(&args, &NoopProgress)?;

    // Verify tracelog was written and contains data
    let csv_content = read_to_string(&tracelog_path)?;
    let lines: Vec<&str> = csv_content.lines().collect();
    assert!(lines.len() >= 2, "Tracelog must have header + at least 1 data row");

    // Verify header contains expected columns
    let header = lines[0];
    assert!(header.contains("n_diff"), "Tracelog header missing n_diff column");
    assert!(header.contains("lh_seq"), "Tracelog header missing lh_seq column");
    assert!(header.contains("lh_pos"), "Tracelog header missing lh_pos column");
    assert!(header.contains("lh_total"), "Tracelog header missing lh_total column");

    // Verify first data row has non-empty likelihood values
    let data_row = lines[1];
    let fields: Vec<&str> = data_row.split(',').collect();
    let columns: Vec<&str> = header.split(',').collect();
    assert!(fields.len() >= 6, "Data row must have at least 6 columns");

    let col = |name: &str| {
      columns
        .iter()
        .position(|c| *c == name)
        .unwrap_or_else(|| panic!("Column '{name}' not found in header"))
    };
    let lh_seq: f64 = fields[col("lh_seq")].parse().expect("lh_seq must be a valid number");
    let lh_pos: f64 = fields[col("lh_pos")].parse().expect("lh_pos must be a valid number");
    let lh_total: f64 = fields[col("lh_total")]
      .parse()
      .expect("lh_total must be a valid number");
    assert!(lh_seq < 0.0, "Sequence log-likelihood must be negative, got {lh_seq}");
    assert!(lh_pos < 0.0, "Positional log-likelihood must be negative, got {lh_pos}");
    assert!(lh_total < 0.0, "Total log-likelihood must be negative, got {lh_total}");

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
