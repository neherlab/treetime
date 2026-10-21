#[cfg(test)]
mod tests {
  use crate::test_utils::{arb_simulated_case, timetree_cell_params};
  use crate::timetree::params::TimetreeParams;
  use helpers::{
    as_test_failure, clade_key_set, internal_nodes_without_valid_interval, run, undated_nodes, younger_parents,
  };
  use proptest::prelude::{ProptestConfig, prop_assert_eq, proptest};
  use treetime_graph::clades::TipIds;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(8))]

    #[test]
    fn test_prop_timetree_simulated_invariants(spec in arb_simulated_case()) {
      let case = spec.simulate().map_err(|error| as_test_failure(&error))?;
      let params = TimetreeParams {
        clock_rate: Some(spec.mu),
        coalescent_opt: true,
        ..timetree_cell_params()
      };

      let output = run(&params, spec.seed, case.timetree_cell_input(spec.length).map_err(|error| as_test_failure(&error))?)?;

      prop_assert_eq!(Vec::<String>::new(), undated_nodes(&output));
      prop_assert_eq!(Vec::<String>::new(), younger_parents(&output));
      let mut ids = TipIds::default();
      let true_keys = clade_key_set(&case.tree.graph, &case.tree.names, &mut ids).map_err(|error| as_test_failure(&error))?;
      let keys = clade_key_set(&output.graph, &output.names, &mut ids).map_err(|error| as_test_failure(&error))?;
      prop_assert_eq!(true_keys, keys);
      prop_assert_eq!(Vec::<String>::new(), internal_nodes_without_valid_interval(&output));
    }
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::error::OperationError;
    use crate::progress::NoopProgress;
    use crate::timetree::params::TimetreeParams;
    use crate::timetree::pipeline::{self, TimetreeInput, TimetreeOutput};
    use eyre::Report;
    use itertools::Itertools;
    use proptest::test_runner::TestCaseError;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::clades::{TipIds, graph_clade_keys};
    use treetime_graph::graph::Graph;
    use treetime_graph::node_names::NodeNames;

    pub(super) fn run(
      params: &TimetreeParams,
      seed: u64,
      input: TimetreeInput,
    ) -> Result<TimetreeOutput, TestCaseError> {
      pipeline::run(params, input, None, None, &NoopCancel, &NoopProgress, &NoopProgress)
        .map_err(OperationError::into_report)
        .map_err(|error| TestCaseError::fail(format!("Timetree run of the case with seed {seed} failed: {error:?}")))
    }

    pub(super) fn as_test_failure(error: &Report) -> TestCaseError {
      TestCaseError::fail(format!("{error:?}"))
    }

    pub(super) fn undated_nodes(output: &TimetreeOutput) -> Vec<String> {
      output
        .graph
        .get_nodes()
        .filter(|node| !output.node_dates.get(&node.key()).is_some_and(Option::is_some))
        .map(|node| output.names[node.key()].to_owned())
        .collect_vec()
    }

    pub(super) fn younger_parents(output: &TimetreeOutput) -> Vec<String> {
      output
        .graph
        .get_edges()
        .filter_map(|edge| {
          let parent = output.node_dates[&edge.source()]?;
          let child = output.node_dates[&edge.target()]?;
          (child < parent).then(|| {
            format!(
              "{} at {child} before its parent {} at {parent}",
              &output.names[edge.target()],
              &output.names[edge.source()]
            )
          })
        })
        .collect_vec()
    }

    pub(super) fn clade_key_set<'a>(
      graph: &Graph,
      names: &'a NodeNames,
      ids: &mut TipIds<'a>,
    ) -> Result<BTreeSet<Vec<usize>>, Report> {
      Ok(graph_clade_keys(graph, names, ids)?.into_values().collect())
    }

    pub(super) fn internal_nodes_without_valid_interval(output: &TimetreeOutput) -> Vec<String> {
      let intervals = output
        .confidence_intervals
        .iter()
        .flatten()
        .map(|interval| (interval.key, (interval.lower, interval.upper)))
        .collect::<BTreeMap<_, _>>();
      output
        .graph
        .get_internal_nodes()
        .filter_map(|node| {
          let name = &output.names[node.key()];
          match intervals.get(&node.key()) {
            None => Some(format!("{name} has no interval")),
            Some(&(lower, upper)) if !(lower.is_finite() && upper.is_finite() && upper > lower) => {
              Some(format!("{name} has the interval [{lower}, {upper}]"))
            },
            Some(_) => None,
          }
        })
        .collect_vec()
    }
  }
}
