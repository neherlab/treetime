#[cfg(test)]
mod tests {
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::solver_outcome::{ConvergenceSummary, UnconvergedEdges};
  use crate::progress::{RunWarning, RunWarningKind};
  use eyre::Report;
  use helpers::{RecordingWarnings, edges_into};
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::{o, vec_of_owned};

  const TREE_NEWICK: &str = "((A:0.1,B:0.1)AB:0.1,C:0.1)root;";

  #[test]
  fn test_solver_outcome_reports_unconverged_rounds_of_a_run_as_one_warning() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut pre_reroot = ConvergenceSummary::default();
    pre_reroot.record(
      &UnconvergedEdges {
        iteration_limit: edges_into(&graph, &names, &["AB", "A"]),
        search_bound: edges_into(&graph, &names, &["C"]),
      },
      &graph,
      &names,
    );
    let mut rounds = ConvergenceSummary::default();
    rounds.record(&UnconvergedEdges::default(), &graph, &names);
    rounds.record(
      &UnconvergedEdges {
        iteration_limit: edges_into(&graph, &names, &["A"]),
        search_bound: vec![],
      },
      &graph,
      &names,
    );
    pre_reroot.merge(rounds);
    let log = RecordingWarnings::default();

    pre_reroot.report(BranchOptMethod::Newton, &log);

    let expected = vec![RunWarning {
      kind: RunWarningKind::BranchLengthsNotConverged,
      message: o!(
        "The newton branch-length optimizer did not converge on 3 edges in 2 of 3 optimization rounds: it reached its iteration limit on 2 edges; it stopped at the upper end of its search interval, where the likelihood still increased, on 1 edge. These edges keep the best length the optimizer evaluated, which may not maximize the likelihood"
      ),
      names: vec_of_owned!["A", "AB", "C"],
    }];
    assert_eq!(expected, log.warnings.into_inner());
    Ok(())
  }

  #[test]
  fn test_solver_outcome_reports_search_bound_alone_in_singular() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut summary = ConvergenceSummary::default();
    summary.record(
      &UnconvergedEdges {
        iteration_limit: vec![],
        search_bound: edges_into(&graph, &names, &["B"]),
      },
      &graph,
      &names,
    );
    let log = RecordingWarnings::default();

    summary.report(BranchOptMethod::Brent, &log);

    let expected = vec![RunWarning {
      kind: RunWarningKind::BranchLengthsNotConverged,
      message: o!(
        "The brent branch-length optimizer did not converge on 1 edge in 1 of 1 optimization round: it stopped at the upper end of its search interval, where the likelihood still increased, on 1 edge. These edges keep the best length the optimizer evaluated, which may not maximize the likelihood"
      ),
      names: vec_of_owned!["B"],
    }];
    assert_eq!(expected, log.warnings.into_inner());
    Ok(())
  }

  #[test]
  fn test_solver_outcome_reports_nothing_when_every_edge_converged() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut summary = ConvergenceSummary::default();
    summary.record(&UnconvergedEdges::default(), &graph, &names);
    summary.record(&UnconvergedEdges::default(), &graph, &names);
    let log = RecordingWarnings::default();

    summary.report(BranchOptMethod::Brent, &log);

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
