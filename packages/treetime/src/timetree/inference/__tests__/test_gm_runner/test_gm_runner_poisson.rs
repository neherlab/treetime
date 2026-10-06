#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{OUTPUTS, load_dates_for_dataset};
  use super::super::test_gm_runner_support::support::{create_poisson_branch_distributions, extract_node_times};
  use crate::clock::date_constraints::load_date_constraints;
  use crate::progress::NoopProgress;
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::bad_branches::{bad_leaves, derive_bad_branches};
  use crate::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::timetree::inference::result::BranchLikelihood;
  use crate::timetree::inference::runner::GRID_POINTS;
  use eyre::Report;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::collections::BTreeSet;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::ebola_20("ebola_20")]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // TODO: enable these datasets when their golden-master gaps are fixed: kb/issues/M-timetree-gm-runner-missing-internal-times.md,
  // kb/issues/M-timetree-date-header-hash.md
  // #[case::dengue_20("dengue_20")]       // TODO: missing internal node times, leaf dates not refined
  // #[case::lassa_l_20("lassa_L_20")]     // TODO: missing internal node times, leaf dates not refined
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // TODO: missing internal node times, leaf dates not refined
  // #[case::rsv_a_20("rsv_a_20")]         // TODO: missing internal node times, leaf dates not refined
  // #[case::tb_20("tb_20")]               // TODO: missing internal node times, leaf dates not refined
  // #[case::zika_20("zika_20")]           // TODO: read_dates strips # from headers, name_column="#name" mismatches
  #[trace]
  #[ignore = "dense-vs-v0 discrepancy: max 0.27 years (grid-width limited, kb/issues/M-timetree-branch-grid-uniform-resolution.md)"]
  fn test_gm_runner_poisson(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = case.poisson();

    let nwk_parsed = nwk_read(case.rerooted_tree_nwk().as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

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
    let bad_branches = derive_bad_branches(&graph, &constraints, &bad_leaves(&graph, &constraints, &BTreeSet::new()))?;
    let backward = propagate_distributions_backward(&graph, &constraints, None, &bad_branches, &branches, MaxGridPoints::default())?;
    let posterior = propagate_distributions_forward(&graph, &constraints, &names, &branches, &backward, MaxGridPoints::default(), &NoopProgress)?;

    let actual = extract_node_times(&graph, &names, &posterior);
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-6);

    Ok(())
  }
}
