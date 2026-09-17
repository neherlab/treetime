#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::cancel::NoopCancel;
  use crate::clock::find_best_root::params::{RerootMethod, RerootSpec};
  use crate::gtr::get_gtr::GtrModelName;
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
  use crate::optimize::pipeline::{OptimizeInput, OptimizeParams, run};
  use crate::progress::NoopProgress;
  use eyre::Report;
  use std::collections::BTreeMap;
  use std::path::Path;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::read_many_fasta_path;
  use treetime_io::nwk::nwk_read_file;
  use treetime_primitives::AlignmentRecord;

  fn load() -> Result<
    (
      Graph,
      BTreeMap<GraphNodeKey, Option<String>>,
      Alphabet,
      Vec<AlignmentRecord>,
      BTreeMap<GraphEdgeKey, Option<f64>>,
    ),
    Report,
  > {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(Path::parent)
      .expect("workspace root");
    let alphabet = Alphabet::default();
    let nwk_parsed = nwk_read_file(workspace_root.join("data/flu/h3n2/20/tree.nwk"))?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln = workspace_root.join("data/flu/h3n2/20/aln.fasta.xz");
    let sequences: Vec<AlignmentRecord> = read_many_fasta_path(&[aln.to_str().expect("utf-8 path")], &alphabet)?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
    Ok((graph, names, alphabet, sequences, branch_lengths))
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
      topology_ops: TopologyOps::default(),
    }
  }

  fn assert_branch_lengths_valid(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) {
    for edge in graph.get_edges() {
      let bl = branch_lengths[&edge.key()].expect("every edge has a branch length after optimization");
      assert!(bl.is_finite() && bl >= 0.0, "invalid branch length {bl}");
    }
  }

  fn root_key(graph: &Graph) -> GraphNodeKey {
    graph.get_exactly_one_root().unwrap().key()
  }

  fn root_child_keys(graph: &Graph) -> Vec<GraphNodeKey> {
    let root = graph.get_exactly_one_root().unwrap();
    root
      .outbound()
      .iter()
      .map(|&edge_key| graph.get_edge(edge_key).unwrap().target())
      .collect()
  }

  fn leaf_names(graph: &Graph, names: &BTreeMap<GraphNodeKey, Option<String>>) -> Vec<String> {
    graph
      .get_leaves()
      .map(|leaf| names.get(&leaf.key()).cloned().flatten().expect("leaf has a name"))
      .collect()
  }

  #[test]
  fn test_optimize_pipeline_reroot_min_dev_changes_root() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let leaves_before = graph.get_leaves().count();
    let root_children_before = root_child_keys(&graph);
    let output = run(
      &params_with(Some(RerootSpec::Method(RerootMethod::MinDev))),
      OptimizeInput {
        graph,
        alphabet,
        sequences,
        branch_lengths,
      },
      &names,
      &NoopCancel,
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().count(), leaves_before);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    let root_children_after = root_child_keys(&output.graph);
    assert_ne!(
      root_children_before, root_children_after,
      "min-dev reroot should change the root position on this unbalanced dataset"
    );
    Ok(())
  }

  #[test]
  fn test_optimize_pipeline_reroot_tips_changes_root() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let leaves_before = graph.get_leaves().count();
    let root_before = root_key(&graph);
    let tips: Vec<String> = leaf_names(&graph, &names).into_iter().take(2).collect();
    let output = run(
      &params_with(Some(RerootSpec::Tips(tips))),
      OptimizeInput {
        graph,
        alphabet,
        sequences,
        branch_lengths,
      },
      &names,
      &NoopCancel,
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().count(), leaves_before);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    let root_after = root_key(&output.graph);
    assert_ne!(root_before, root_after, "tip-based reroot should move the root");
    Ok(())
  }

  #[test]
  fn test_optimize_pipeline_reroot_min_dev_dense_completes() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let leaves_before = graph.get_leaves().count();

    let mut params = params_with(Some(RerootSpec::Method(RerootMethod::MinDev)));
    params.dense = Some(true);
    let output = run(
      &params,
      OptimizeInput {
        graph,
        alphabet,
        sequences,
        branch_lengths,
      },
      &names,
      &NoopCancel,
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().count(), leaves_before);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    Ok(())
  }

  // Keeping the root (the default) leaves the run unchanged relative to no reroot
  // flags; both must complete and preserve the tree.
  #[test]
  fn test_optimize_pipeline_keep_root_completes() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let leaves_before = graph.get_leaves().count();
    let output = run(
      &params_with(None),
      OptimizeInput {
        graph,
        alphabet,
        sequences,
        branch_lengths,
      },
      &names,
      &NoopCancel,
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().count(), leaves_before);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    Ok(())
  }
}
