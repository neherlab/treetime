#[cfg(test)]
mod tests {
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output::{LadderizeArg, OutputCoreArgs, TimetreeOutputSelection, TopologyOrderArgs};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::output::coalescent::{CoalescentOutput, CoalescentOutputMode};
  use crate::commands::timetree::run::run_timetree_estimation;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use std::fs::read_to_string;
  use std::path::PathBuf;
  use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
  use treetime_utils::io::json::json_read_file;

  #[test]
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
  fn test_pipeline_timetree_convergence() -> Result<(), Report> {
    let root = project_root();
    let outdir = root.join("tmp/test-convergence-pipeline");
    std::fs::create_dir_all(&outdir)?;
    let tracelog_path = outdir.join("tracelog.csv");

    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      max_iter: 3,
      output: OutputCoreArgs {
        output_all: Some(outdir.clone()),
        ..Default::default()
      },
      output_tracelog: Some(tracelog_path.clone()),
      ..TreetimeTimetreeArgsRaw::default()
    })
    .unwrap();

    run_timetree_estimation(&args, &NoopProgress)?;

    // Verify tracelog was written and contains data
    let csv_content = read_to_string(&tracelog_path)?;
    let lines: Vec<&str> = csv_content.lines().collect();
    assert!(lines.len() >= 2, "Tracelog must have header + at least 1 data row");

    // Verify the exact log-likelihood schema.
    let header = lines[0];
    let expected_header =
      "n_diff,n_resolved,max_time_change,rms_time_change,log_lh_seq,log_lh_pos,log_lh_coal,log_lh_total";
    assert_eq!(expected_header, header);

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
    let log_lh_seq: f64 = fields[col("log_lh_seq")]
      .parse()
      .expect("log_lh_seq must be a valid number");
    let log_lh_pos: f64 = fields[col("log_lh_pos")]
      .parse()
      .expect("log_lh_pos must be a valid number");
    let log_lh_total: f64 = fields[col("log_lh_total")]
      .parse()
      .expect("log_lh_total must be a valid number");
    assert!(
      log_lh_seq < 0.0,
      "Sequence log-likelihood must be negative, got {log_lh_seq}"
    );
    assert!(
      log_lh_pos < 0.0,
      "Positional log-likelihood must be negative, got {log_lh_pos}"
    );
    assert!(
      log_lh_total < 0.0,
      "Total log-likelihood must be negative, got {log_lh_total}"
    );

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

  #[test]
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
  fn test_pipeline_timetree_ladderize_applies_to_auspice() -> Result<(), Report> {
    let root = project_root();
    let output = tempfile::tempdir()?;
    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      max_iter: 0,
      output: OutputCoreArgs {
        output_all: Some(output.path().to_path_buf()),
        ..OutputCoreArgs::default()
      },
      output_selection: vec![TimetreeOutputSelection::Auspice],
      topology_order: TopologyOrderArgs {
        ladderize: Some(LadderizeArg::Descending),
        ..TopologyOrderArgs::default()
      },
      ..TreetimeTimetreeArgsRaw::default()
    })
    .unwrap();

    run_timetree_estimation(&args, &NoopProgress)?;

    let tree: AuspiceTree = json_read_file(output.path().join("timetree.auspice.json"))?;
    let actual = tree.tree.children.iter().map(count_leaves).collect::<Vec<_>>();
    let mut expected = actual.clone();
    expected.sort_unstable_by(|left, right| right.cmp(left));
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  #[ignore = "mass-sized node times break downstream invariants (positional log-lh, polytomy resolution): kb/issues/H-timetree-mass-sizing-node-times-break-downstream-invariants.md"]
  fn test_pipeline_timetree_writes_coalescent_outputs() -> Result<(), Report> {
    let root = project_root();
    let output = tempfile::tempdir()?;
    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      max_iter: 3,
      coalescent_skyline: true,
      output: OutputCoreArgs {
        output_all: Some(output.path().to_path_buf()),
        ..Default::default()
      },
      output_coalescent_json: Some(output.path().join("timetree.coalescent.json")),
      ..TreetimeTimetreeArgsRaw::default()
    })
    .unwrap();

    run_timetree_estimation(&args, &NoopProgress)?;

    // The default `--output-all` set writes the coalescent TSV; the JSON is requested per file.
    let tsv_path = output.path().join("timetree.coalescent.tsv");
    assert!(tsv_path.exists(), "default --output-all must write the coalescent TSV");

    let doc: CoalescentOutput = json_read_file(output.path().join("timetree.coalescent.json"))?;
    assert_eq!(CoalescentOutputMode::Skyline, doc.inputs.mode);
    assert!(
      !doc.outputs.segments.is_empty(),
      "a skyline run must report at least one coalescent segment"
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

  fn count_leaves(node: &AuspiceTreeNode) -> usize {
    if node.children.is_empty() {
      1
    } else {
      node.children.iter().map(count_leaves).sum()
    }
  }
}
