#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::project_root;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::metadata::MetadataArgs;
  use crate::commands::shared::output_args::TimetreeOutputSelectionArgs;
  use crate::commands::shared::output_args::{OutputCoreArgs, TimetreeOutputSelection};
  use crate::commands::shared::seed::SeedArgs;
  use crate::commands::shared::topology_order_args::{TopologyOrderArg, TopologyOrderArgs};
  use crate::commands::shared::tree_input::TreeArgs;
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::coalescent::{CoalescentArg, parse_coalescent};
  use crate::commands::timetree::run::run_timetree_estimation;
  use eyre::{Report, WrapErr};
  use helpers::{coalescent_node_dates, count_leaves};
  use std::collections::BTreeMap;
  use std::fs::read_to_string;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::auspice_types::{AuspiceTree, AuspiceTreeNode};
  use treetime_utils::io::json::{json_read_file, json_value_read_file};
  use treetime_utils::make_report;

  #[test]
  fn test_pipeline_timetree_convergence() -> Result<(), Report> {
    let root = project_root();
    let outdir = root.join("tmp/test-convergence-pipeline");
    std::fs::create_dir_all(&outdir).wrap_err("When creating the output directory")?;
    let tracelog_path = outdir.join("tracelog.csv");

    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: TreeArgs {
        tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      },
      metadata: MetadataArgs {
        metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      },
      max_iter: 3,
      output: OutputCoreArgs {
        output_all: Some(outdir.clone()),
        ..Default::default()
      },
      output_tracelog: Some(tracelog_path.clone()),
      seed_args: SeedArgs { seed: Some(7) },
      ..TreetimeTimetreeArgsRaw::default()
    })
    .unwrap();

    run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let csv_content = read_to_string(&tracelog_path).wrap_err("When reading the trace log")?;
    let lines: Vec<&str> = csv_content.lines().collect();
    assert!(lines.len() >= 2, "Tracelog must have header + at least 1 data row");

    let header = lines[0];
    let expected_header = "n_diff,n_resolved,max_time_change,rms_time_change,log_lh_seq,log_lh_time,log_lh_coal";
    assert_eq!(expected_header, header);

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
    let log_lh_time: f64 = fields[col("log_lh_time")]
      .parse()
      .expect("log_lh_time must be a valid number");
    assert!(
      log_lh_seq < 0.0,
      "Sequence log-likelihood must be negative, got {log_lh_seq}"
    );
    assert!(
      log_lh_time.is_finite(),
      "Node-time log-likelihood must be finite, got {log_lh_time}"
    );

    let nwk_path = outdir.join("timetree.nwk");
    let nex_path = outdir.join("timetree.nexus");
    assert!(nwk_path.exists(), "Output newick file must exist");
    assert!(nex_path.exists(), "Output nexus file must exist");
    assert!(
      std::fs::metadata(&nwk_path)
        .wrap_err("When reading the Newick output metadata")?
        .len()
        > 0,
      "Output newick file must be non-empty"
    );

    Ok(())
  }

  #[test]
  fn test_pipeline_timetree_topology_order_applies_to_auspice() -> Result<(), Report> {
    let root = project_root();
    let output = tempfile::tempdir().wrap_err("When creating a temporary directory")?;
    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: TreeArgs {
        tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      },
      metadata: MetadataArgs {
        metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      },
      max_iter: 0,
      output: OutputCoreArgs {
        output_all: Some(output.path().to_path_buf()),
        ..OutputCoreArgs::default()
      },
      output_selection: TimetreeOutputSelectionArgs {
        output_selection: vec![TimetreeOutputSelection::Auspice],
      },
      topology_order: TopologyOrderArgs {
        topology_order: Some(TopologyOrderArg::DescendantCountReverse),
        ..TopologyOrderArgs::default()
      },
      seed_args: SeedArgs { seed: Some(7) },
      ..TreetimeTimetreeArgsRaw::default()
    })
    .unwrap();

    run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let tree: AuspiceTree = json_read_file(output.path().join("timetree.auspice.json"))?;
    let actual = tree.tree.children.iter().map(count_leaves).collect::<Vec<_>>();
    let mut expected = actual.clone();
    expected.sort_unstable_by(|left, right| right.cmp(left));
    assert_eq!(expected, actual);

    Ok(())
  }

  #[test]
  fn test_pipeline_timetree_writes_coalescent_outputs() -> Result<(), Report> {
    let root = project_root();
    let output = tempfile::tempdir().wrap_err("When creating a temporary directory")?;
    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: TreeArgs {
        tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
      },
      metadata: MetadataArgs {
        metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
      },
      max_iter: 3,
      coalescent: Some(CoalescentArg::Skyline),
      output: OutputCoreArgs {
        output_all: Some(output.path().to_path_buf()),
        ..Default::default()
      },
      output_coalescent_json: Some(output.path().join("timetree.coalescent.json")),
      seed_args: SeedArgs { seed: Some(7) },
      ..TreetimeTimetreeArgsRaw::default()
    })
    .unwrap();

    run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let tsv_path = output.path().join("timetree.coalescent.tsv");
    assert!(tsv_path.exists(), "default --output-all must write the coalescent TSV");

    let doc = json_value_read_file(output.path().join("timetree.coalescent.json"))?;
    assert_eq!("skyline", doc["inputs"]["mode"]);
    assert!(
      doc["outputs"]["segments"]
        .as_array()
        .is_some_and(|segments| !segments.is_empty()),
      "a skyline run must report at least one coalescent segment"
    );

    Ok(())
  }

  #[test]
  fn test_pipeline_timetree_coalescent_in_substitutions_per_site_equals_its_time_at_a_fixed_clock_rate()
  -> Result<(), Report> {
    let clock_rate = 0.004;
    let substitutions = 0.01;

    let in_substitutions = coalescent_node_dates(&substitutions.to_string(), Some(clock_rate))?;
    let in_years = coalescent_node_dates(&format!("{}y", substitutions / clock_rate), Some(clock_rate))?;

    assert_eq!(in_years, in_substitutions);
    Ok(())
  }

  #[test]
  fn test_pipeline_timetree_coalescent_in_years_keeps_the_dates_of_a_fixed_time_scale() -> Result<(), Report> {
    let expected: BTreeMap<String, f64> = json_read_file(
      project_root()
        .join("packages/app-commands/src/commands/timetree/__tests__/__fixtures__/coalescent_years_node_dates.json"),
    )?;

    let mut actual = coalescent_node_dates("2.56y", None)?;
    actual.retain(|name, _| expected.contains_key(name));

    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn count_leaves(node: &AuspiceTreeNode) -> usize {
      if node.children.is_empty() {
        1
      } else {
        node.children.iter().map(count_leaves).sum()
      }
    }

    pub(super) fn coalescent_node_dates(
      coalescent: &str,
      clock_rate: Option<f64>,
    ) -> Result<BTreeMap<String, f64>, Report> {
      let root = project_root();
      let output = tempfile::tempdir().wrap_err("When creating a temporary directory")?;
      let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
        },
        tree: TreeArgs {
          tree: Some(root.join("data/flu/h3n2/20/tree.nwk")),
        },
        metadata: MetadataArgs {
          metadata: Some(root.join("data/flu/h3n2/20/metadata.tsv")),
        },
        max_iter: 2,
        clock_rate,
        coalescent: Some(parse_coalescent(coalescent)?),
        output: OutputCoreArgs {
          output_all: Some(output.path().to_path_buf()),
          ..Default::default()
        },
        seed_args: SeedArgs { seed: Some(7) },
        ..TreetimeTimetreeArgsRaw::default()
      })?;

      run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

      let node_data = json_value_read_file(output.path().join("timetree.node-data.json"))?;
      let nodes = node_data["nodes"]
        .as_object()
        .ok_or_else(|| make_report!("The augur node data has no nodes"))?;
      nodes
        .iter()
        .map(|(name, node)| {
          let date = node["numdate"]
            .as_f64()
            .ok_or_else(|| make_report!("Node {name} has no numdate"))?;
          Ok((name.clone(), date))
        })
        .collect()
    }
  }
}
