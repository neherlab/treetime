#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use crate::gtr::get_gtr::GtrModelName;
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
  use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
  use crate::payload::ancestral::GraphAncestral;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use std::path::Path;
  use treetime_graph::edge::HasBranchLength;
  use treetime_graph::node::{GraphNodeKey, Named};
  use treetime_io::fasta::{FastaRecord, read_many_fasta};
  use treetime_io::nwk::nwk_read_file;

  fn load() -> Result<(GraphAncestral, Alphabet, Vec<FastaRecord>), Report> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(Path::parent)
      .expect("workspace root");
    let alphabet = Alphabet::default();
    let graph: GraphAncestral = nwk_read_file(workspace_root.join("data/flu/h3n2/20/tree.nwk"))?;
    let aln = workspace_root.join("data/flu/h3n2/20/aln.fasta.xz");
    let sequences = read_many_fasta(&[aln.to_str().expect("utf-8 path")], &alphabet)?;
    Ok((graph, alphabet, sequences))
  }

  fn params_with(reroot_spec: Option<RerootSpec>) -> OptimizeParams {
    OptimizeParams {
      model: GtrModelName::JC69,
      dense: Some(false),
      max_iter: 2,
      dp: 0.1,
      damping: 0.75,
      opt_method: BranchOptMethod::default(),
      initial_guess: InitialGuessMode::default(),
      no_indels: false,
      reroot_spec,
    }
  }

  fn assert_branch_lengths_valid(graph: &GraphAncestral) {
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let bl = edge
        .payload()
        .read_arc()
        .branch_length()
        .expect("every edge has a branch length after optimization");
      assert!(bl.is_finite() && bl >= 0.0, "invalid branch length {bl}");
    }
  }

  fn root_key(graph: &GraphAncestral) -> GraphNodeKey {
    graph.get_exactly_one_root().unwrap().read_arc().key()
  }

  fn root_child_keys(graph: &GraphAncestral) -> Vec<GraphNodeKey> {
    let root = graph.get_exactly_one_root().unwrap();
    let root = root.read_arc();
    root
      .outbound()
      .iter()
      .map(|&edge_key| graph.get_edge(edge_key).unwrap().read_arc().target())
      .collect()
  }

  fn leaf_names(graph: &GraphAncestral) -> Vec<String> {
    graph
      .get_leaves()
      .iter()
      .map(|leaf| {
        leaf
          .read_arc()
          .payload()
          .read_arc()
          .name()
          .map(|name| name.as_ref().to_owned())
          .expect("leaf has a name")
      })
      .collect()
  }

  #[test]
  fn test_optimize_pipeline_reroot_min_dev_changes_root() -> Result<(), Report> {
    let (graph, alphabet, sequences) = load()?;
    let leaves_before = graph.get_leaves().len();
    let root_children_before = root_child_keys(&graph);

    let output = run(
      &params_with(Some(RerootSpec::Method(RerootMethod::MinDev))),
      OptimizeInput {
        graph,
        alphabet,
        sequences,
      },
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().len(), leaves_before);
    assert_branch_lengths_valid(&output.graph);
    let root_children_after = root_child_keys(&output.graph);
    assert_ne!(
      root_children_before, root_children_after,
      "min-dev reroot should change the root position on this unbalanced dataset"
    );
    Ok(())
  }

  #[test]
  fn test_optimize_pipeline_reroot_tips_changes_root() -> Result<(), Report> {
    let (graph, alphabet, sequences) = load()?;
    let leaves_before = graph.get_leaves().len();
    let root_before = root_key(&graph);
    let tips: Vec<String> = leaf_names(&graph).into_iter().take(2).collect();

    let output = run(
      &params_with(Some(RerootSpec::Tips(tips))),
      OptimizeInput {
        graph,
        alphabet,
        sequences,
      },
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().len(), leaves_before);
    assert_branch_lengths_valid(&output.graph);
    let root_after = root_key(&output.graph);
    assert_ne!(root_before, root_after, "tip-based reroot should move the root");
    Ok(())
  }

  #[test]
  fn test_optimize_pipeline_reroot_min_dev_dense_completes() -> Result<(), Report> {
    let (graph, alphabet, sequences) = load()?;
    let leaves_before = graph.get_leaves().len();

    let mut params = params_with(Some(RerootSpec::Method(RerootMethod::MinDev)));
    params.dense = Some(true);

    let output = run(
      &params,
      OptimizeInput {
        graph,
        alphabet,
        sequences,
      },
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().len(), leaves_before);
    assert_branch_lengths_valid(&output.graph);
    Ok(())
  }

  // Keeping the root (the default) leaves the run unchanged relative to no reroot
  // flags; both must complete and preserve the tree.
  #[test]
  fn test_optimize_pipeline_keep_root_completes() -> Result<(), Report> {
    let (graph, alphabet, sequences) = load()?;
    let leaves_before = graph.get_leaves().len();

    let output = run(
      &params_with(None),
      OptimizeInput {
        graph,
        alphabet,
        sequences,
      },
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().len(), leaves_before);
    assert_branch_lengths_valid(&output.graph);
    Ok(())
  }
}
