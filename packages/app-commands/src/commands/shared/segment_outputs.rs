use app_output::annotated_graph::AnnotatedGraph;
use app_output::output_plan::{OutputSelection, ResolvedOutputs, TreeWriteKind};
use app_output::tree_output::write_graph_outputs;
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;
use treetime::progress::LogSink;
use treetime::{progress_info, progress_warn};
use treetime_graph::graph::Graph;
use treetime_graph::scope::{EdgeScopes, Inheritance};

pub(crate) fn segment_graphs(
  graph: &Graph,
  scopes: Option<&Arc<EdgeScopes>>,
) -> Result<Option<BTreeMap<usize, Graph>>, Report> {
  let Some(scopes) = scopes else {
    return Ok(None);
  };
  (0..scopes.segment_count())
    .map(|segment| {
      let inheritance = Inheritance::segment(Arc::clone(scopes), segment)?;
      Ok((segment + 1, inheritance.subgraph(graph)?))
    })
    .collect::<Result<BTreeMap<_, _>, Report>>()
    .map(Some)
}

pub(crate) fn write_whole_graph_outputs(
  graph: &AnnotatedGraph<'_>,
  resolved: &ResolvedOutputs,
  run_files: &[OutputSelection],
  log: &dyn LogSink,
) -> Result<(), Report> {
  let (written, skipped): (Vec<_>, Vec<_>) = resolved
    .tree_outputs
    .iter()
    .map(|(&kind, path)| (kind, path.clone()))
    .partition(|(kind, _)| matches!(kind, TreeWriteKind::GraphJson | TreeWriteKind::Dot));
  for (kind, path) in written {
    write_graph_outputs(graph, &BTreeMap::from([(kind, path.clone())]))?;
    progress_info!(log, "Wrote graph to {}", path.display());
  }
  let skipped: Vec<PathBuf> = skipped
    .into_iter()
    .map(|(_, path)| path)
    .chain(
      run_files
        .iter()
        .filter_map(|&selection| resolved.path(selection))
        .map(PathBuf::from),
    )
    .collect();
  if !skipped.is_empty() {
    progress_warn!(
      log,
      "Not writing {}: these files hold one tree for the whole run, but --segment-trees gives a graph of several segment trees. Put {{alignment}} into a tree output path to write one file per segment",
      skipped.iter().map(|path| format!("'{}'", path.display())).join(", ")
    );
  }
  Ok(())
}
