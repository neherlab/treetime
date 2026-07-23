#[cfg(test)]
mod tests {
  use crate::coalescent::total_lh::compute_coalescent_total_lh;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output::{LadderizeArg, OutputCoreArgs, TimetreeOutputSelection, TopologyOrderArgs};
  use crate::commands::timetree::args::TreetimeTimetreeArgs;
  use crate::commands::timetree::initialization::load_input_data;
  use crate::commands::timetree::run::{run_timetree_estimation, timetree_params_from_args};
  use crate::progress::NoopProgress;
  use crate::timetree::pipeline::{self, TimetreeInput};
  use eyre::Report;
  use std::fs::read_to_string;
  use std::path::PathBuf;
  use treetime_distribution::Distribution;
  use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
  use treetime_utils::io::json::json_read_file;

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
      output: OutputCoreArgs {
        output_all: Some(outdir.clone()),
        ..Default::default()
      },
      output_tracelog: Some(tracelog_path.clone()),
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

  #[test]
  fn test_pipeline_timetree_ladderize_applies_to_auspice() -> Result<(), Report> {
    let root = project_root();
    let output = tempfile::tempdir()?;
    let args = TreetimeTimetreeArgs {
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
      ..TreetimeTimetreeArgs::default()
    };

    run_timetree_estimation(&args, &NoopProgress)?;

    let tree: AuspiceTree = json_read_file(output.path().join("timetree.auspice.json"))?;
    let actual = tree.tree.children.iter().map(count_leaves).collect::<Vec<_>>();
    let mut expected = actual.clone();
    expected.sort_unstable_by(|left, right| right.cmp(left));
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_pipeline_coalescent_opt_prior_maximizes_output_tree() -> Result<(), Report> {
    // B1 regression: the reported coalescent prior must be fit on the final output
    // node times, so it must maximize the output tree's coalescent likelihood. A
    // prior fit on the first, least-refined tree (the pre-fix behavior) would not.
    let root = project_root();
    let args = TreetimeTimetreeArgs {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      max_iter: 2,
      coalescent_opt: true,
      ..TreetimeTimetreeArgs::default()
    };

    let input_data = load_input_data(&args)?;
    let params = timetree_params_from_args(&args);
    let input = TimetreeInput {
      graph: input_data.graph,
      alphabet: input_data.alphabet,
      sequences: input_data.aln,
      dates: input_data.dates,
    };

    let output = pipeline::run(&params, input, None, &NoopProgress)?;

    let tc_dist = output.coalescent_tc.expect("coalescent-opt must report a fitted Tc");
    // A constant prior evaluates to the same Tc across its domain.
    let tc = tc_dist.eval(2005.0)?;
    assert!(
      tc.is_finite() && tc > 0.0,
      "reported Tc must be finite and positive, got {tc}"
    );

    // Optimality oracle (independent of optimize_tc): the coalescent log-likelihood
    // L(Tc) = -M ln Tc - I/Tc is strictly unimodal, so the MLE of the output tree
    // beats any perturbed Tc. This fails if the prior was fit on a different tree.
    let lh_at = |scale: f64| compute_coalescent_total_lh(&output.graph, &Distribution::constant(tc * scale));
    let lh_opt = lh_at(1.0)?;
    let lh_lo = lh_at(0.8)?;
    let lh_hi = lh_at(1.25)?;
    assert!(
      lh_opt >= lh_lo && lh_opt >= lh_hi,
      "reported Tc={tc} must maximize the output-tree coalescent likelihood: \
       L(Tc)={lh_opt}, L(0.8·Tc)={lh_lo}, L(1.25·Tc)={lh_hi}"
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
