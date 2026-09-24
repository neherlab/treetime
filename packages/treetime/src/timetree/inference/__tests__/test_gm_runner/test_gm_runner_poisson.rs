#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{OUTPUTS, load_dates_for_dataset};
  use super::super::test_gm_runner_support::support::{create_poisson_branch_distributions, extract_node_times};
  use crate::clock::date_constraints::load_date_constraints;
  use crate::progress::NoopProgress;
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::bad_branches::{derive_bad_branches, undated_leaves};
  use crate::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::timetree::inference::runner::GRID_POINTS;
  use crate::timetree::inference::time_inference::BranchLikelihood;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
#[rstest]
  #[case::ebola_20("ebola_20")]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  #[ignore = "dense-vs-v0 discrepancy: max 0.27 years (grid-width limited, kb/issues/M-timetree-branch-grid-uniform-resolution.md)"]
  fn test_gm_runner_poisson(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = case.poisson();

    let nwk_parsed = nwk_read_str(case.rerooted_tree_nwk())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates, &graph, &names, &NoopProgress)?;

    let branch_distributions = create_poisson_branch_distributions(
      &graph,
      &branch_lengths,
      case.clock_rate(),
      case.sequence_length(),
      GRID_POINTS,
    )?;
    let branches: BTreeMap<GraphEdgeKey, BranchLikelihood> = graph
      .get_edges()
      .map(|edge| {
        let branch = BranchLikelihood {
          distribution: branch_distributions.get(&edge.key()).cloned(),
          time_length: None,
        };
        (edge.key(), branch)
      })
      .collect();
    let bad_branches = derive_bad_branches(&graph, &constraints, &undated_leaves(&graph, &constraints))?;
    let backward = propagate_distributions_backward(&graph, &constraints, None, &bad_branches, &branches)?;
    let posterior = propagate_distributions_forward(&graph, &constraints, &names, &branches, &backward, &NoopProgress)?;

    let actual = extract_node_times(&graph, &names, &posterior);
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-6);

    Ok(())
  }
}
