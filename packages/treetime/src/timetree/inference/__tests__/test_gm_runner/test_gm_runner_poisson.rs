#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{OUTPUTS, load_dates_for_dataset};
  use crate::clock::date_constraints::load_date_constraints;
  use crate::partition::timetree::GraphTimetree;
  use crate::timetree::inference::backward_pass::propagate_distributions_backward;
  use crate::timetree::inference::forward_pass::propagate_distributions_forward;
  use crate::timetree::inference::runner::BRANCH_GRID_SIZE;
  use crate::timetree::utils::{create_poisson_branch_distributions, extract_node_times};
  use eyre::Report;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  // --- Poisson tests ---

  #[rustfmt::skip]
#[rstest]
  // #[case::dengue_20("dengue_20")]       // TODO: missing internal node times, leaf dates not refined
  #[case::ebola_20("ebola_20")]
  #[case::flu_h3n2_20("flu_h3n2_20")]
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

    let graph: GraphTimetree = nwk_read_str(case.rerooted_tree_nwk())?;
    let dates = load_dates_for_dataset(dataset)?;
    load_date_constraints(&dates, &graph)?;

    create_poisson_branch_distributions(
      &graph,
      case.clock_rate(),
      case.sequence_length(),
      BRANCH_GRID_SIZE,
    )?;
    propagate_distributions_backward(&graph, None)?;
    propagate_distributions_forward(&graph)?;

    let actual = extract_node_times(&graph);
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-6);

    Ok(())
  }
}
