#[cfg(test)]
mod tests {
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::solver_outcome::UnconvergedEdges;
  use crate::progress::{RunWarning, RunWarningKind};
  use eyre::Report;
  use helpers::{RecordingWarnings, edges_into};
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::{o, vec_of_owned};

  const TREE_NEWICK: &str = "((A:0.1,B:0.1)AB:0.1,C:0.1)root;";

  #[test]
  fn test_solver_outcome_reports_unconverged_edges_as_one_warning() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let unconverged = UnconvergedEdges {
      edges: edges_into(&graph, &names, &["AB", "A"]),
      optimized: 4,
    };
    let log = RecordingWarnings::default();

    unconverged.report(BranchOptMethod::Newton, &graph, &names, &log);

    let expected = vec![RunWarning {
      kind: RunWarningKind::BranchLengthsNotConverged,
      message: o!(
        "The newton branch-length optimizer reached its iteration limit on 2 of 4 edges. These edges keep the best length the optimizer evaluated, which may not maximize the likelihood"
      ),
      names: vec_of_owned!["A", "AB"],
    }];
    assert_eq!(expected, log.warnings.into_inner());
    Ok(())
  }

  #[test]
  fn test_solver_outcome_reports_nothing_when_every_edge_converged() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let unconverged = UnconvergedEdges {
      edges: vec![],
      optimized: 4,
    };
    let log = RecordingWarnings::default();

    unconverged.report(BranchOptMethod::Brent, &graph, &names, &log);

    assert_eq!(Vec::<RunWarning>::new(), log.warnings.into_inner());
    Ok(())
  }

  mod helpers {
    use crate::progress::{LogLevel, LogSink, RunWarning};
    use parking_lot::Mutex;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;

    #[derive(Default)]
    pub(super) struct RecordingWarnings {
      pub(super) warnings: Mutex<Vec<RunWarning>>,
    }

    impl LogSink for RecordingWarnings {
      fn log(&self, _level: LogLevel, _message: &str) {}

      fn log_enabled(&self, _level: LogLevel) -> bool {
        true
      }

      fn warning(&self, warning: &RunWarning) {
        self.warnings.lock().push(warning.clone());
      }
    }

    pub(super) fn edges_into(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      targets: &[&str],
    ) -> Vec<GraphEdgeKey> {
      targets
        .iter()
        .map(|&target| {
          graph
            .get_edges()
            .find(|edge| names[&edge.target()].as_deref() == Some(target))
            .expect("edge into the named node")
            .key()
        })
        .collect()
    }
  }
}
