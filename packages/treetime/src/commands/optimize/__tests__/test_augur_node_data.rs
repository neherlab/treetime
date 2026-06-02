#[cfg(test)]
mod tests {
  use approx::assert_relative_eq;
  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::json_read_str;
  use util_augur_node_data_json::AugurNodeDataJsonRefine;

  // optimize node data is the augur non-timetree refine shape: per-node
  // branch_length (the ML-optimized divergence, subs/site) plus alignment and
  // input_tree metadata. No clock, date, or confidence fields.

  #[test]
  fn test_augur_node_data_optimize_branch_lengths() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");

    // TODO(kb/issues/N-nwk-branch-length-f32-precision-loss.md): the `bio` crate
    // Newick reader parses weights as f32; leaf values match only to f32 relative
    // precision after widening to f64. Root is set to exactly 0.0 by the writer.
    assert_relative_eq!(data.nodes["leaf_a"].branch_length, 0.005, max_relative = 1e-6);
    assert_relative_eq!(data.nodes["leaf_b"].branch_length, 0.010, max_relative = 1e-6);
    assert_relative_eq!(data.nodes["root"].branch_length, 0.0);
  }

  #[test]
  fn test_augur_node_data_optimize_only_branch_length_field() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");

    // Every timetree-only and confidence field must be omitted: optimize does no
    // temporal inference and v1 does not parse input-tree branch support.
    for (name, node) in &data.nodes {
      assert!(node.confidence.is_none(), "{name}: confidence must be omitted");
      assert!(node.numdate.is_none(), "{name}: numdate must be omitted");
      assert!(node.clock_length.is_none(), "{name}: clock_length must be omitted");
      assert!(
        node.mutation_length.is_none(),
        "{name}: mutation_length must be omitted"
      );
      assert!(node.raw_date.is_none(), "{name}: raw_date must be omitted");
      assert!(node.date.is_none(), "{name}: date must be omitted");
      assert!(node.date_inferred.is_none(), "{name}: date_inferred must be omitted");
      assert!(
        node.num_date_confidence.is_none(),
        "{name}: num_date_confidence must be omitted"
      );
    }
  }

  #[test]
  fn test_augur_node_data_optimize_metadata() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");

    // clock is absent (no time-tree inference); alignment and input_tree carry
    // the input file paths passed to the writer.
    assert!(data.metadata.clock.is_none());
    assert_eq!(data.metadata.alignment.as_deref(), Some("aln.fasta"));
    assert_eq!(data.metadata.input_tree.as_deref(), Some("tree.nwk"));
  }

  #[test]
  fn test_augur_node_data_optimize_generated_by() {
    let data = helpers::write_and_read("(leaf_a:0.005,leaf_b:0.010)root;");
    let generated_by = data.generated_by.unwrap();
    assert_eq!(generated_by.program, "treetime");
    assert_eq!(generated_by.version, env!("CARGO_PKG_VERSION"));
  }

  #[test]
  fn test_augur_node_data_optimize_roundtrip() {
    let json_str = helpers::write_json("(leaf_a:0.005,leaf_b:0.010)root;");

    let original: serde_json::Value = serde_json::from_str(&json_str).unwrap();
    let typed: AugurNodeDataJsonRefine = json_read_str(&json_str).unwrap();
    let roundtripped: serde_json::Value = serde_json::to_value(&typed).unwrap();

    assert_eq!(original, roundtripped);
  }

  // End-to-end: the optimize command writes a readable node data file whose
  // branch lengths are finite and non-negative after optimization.
  #[test]
  fn test_augur_node_data_optimize_end_to_end() {
    use crate::commands::optimize::args::TreetimeOptimizeArgs;
    use crate::commands::optimize::run::run_optimize;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::shared::output::OutputArgs;
    use crate::progress::NoopProgress;

    let root = helpers::project_root();
    let outdir = root.join("tmp/test-optimize-node-data");
    std::fs::create_dir_all(&outdir).unwrap();

    let args = TreetimeOptimizeArgs {
      alignment: AlignmentArgs {
        alignment: vec![root.join("data/flu/h3n2/20/aln.fasta.xz")],
      },
      tree: root.join("data/flu/h3n2/20/tree.nwk"),
      max_iter: 2,
      output: OutputArgs {
        outdir: outdir.clone(),
        ..Default::default()
      },
      ..TreetimeOptimizeArgs::default()
    };

    run_optimize(&args, &NoopProgress).unwrap();

    let json_str = std::fs::read_to_string(outdir.join("optimize.augur-node-data.json")).unwrap();
    let data: AugurNodeDataJsonRefine = json_read_str(&json_str).unwrap();

    assert!(!data.nodes.is_empty(), "node data must contain nodes");
    for (name, node) in &data.nodes {
      assert!(
        node.branch_length.is_finite() && node.branch_length >= 0.0,
        "{name}: branch_length must be finite and non-negative, got {}",
        node.branch_length
      );
    }
  }

  mod helpers {
    use crate::commands::optimize::augur_node_data::write_augur_node_data_json;
    use crate::payload::ancestral::GraphAncestral;
    use std::path::{Path, PathBuf};
    use tempfile::NamedTempFile;
    use treetime_io::nwk::nwk_read_str;
    use treetime_utils::io::json::json_read_str;
    use util_augur_node_data_json::AugurNodeDataJsonRefine;

    pub fn write_json(nwk: &str) -> String {
      let graph: GraphAncestral = nwk_read_str(nwk).unwrap();
      let tmp = NamedTempFile::new().unwrap();
      write_augur_node_data_json(
        &graph,
        Some(Path::new("aln.fasta")),
        Some(Path::new("tree.nwk")),
        tmp.path(),
      )
      .unwrap();
      std::fs::read_to_string(tmp.path()).unwrap()
    }

    pub fn write_and_read(nwk: &str) -> AugurNodeDataJsonRefine {
      json_read_str(write_json(nwk)).unwrap()
    }

    pub fn project_root() -> PathBuf {
      PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .map(PathBuf::from)
        .expect("project has workspace root")
    }
  }
}
