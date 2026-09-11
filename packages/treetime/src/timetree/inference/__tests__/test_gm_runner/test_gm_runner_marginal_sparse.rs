#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };

  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::initialize_marginal;
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_state::ClockState;
  use crate::clock::date_constraints::load_date_constraints;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::timetree::inference::runner::run_timetree;
  use crate::timetree::timetree_state::TimetreeState;
  use crate::timetree::utils::{extract_node_times, initialize_node_divergences};

  use crate::partition::timetree::partition::{GraphTimetree, PartitionTimetree};
  use eyre::Report;

  use rstest::rstest;

  use treetime_io::nwk::{NwkParse, nwk_read_str};

  use treetime_utils::pretty_assert_map_abs_diff_eq;

  // Cross-mode consistency: sparse marginal inference should produce node times
  // matching dense marginal inference. No independent sparse oracle exists (v0 has
  // no sparse mode), so dense golden master output serves as reference.

  #[rustfmt::skip]
#[rstest]
  // #[case::dengue_20("dengue_20")]       // TODO: missing internal node times, leaf dates not refined
  // #[case::ebola_20("ebola_20")]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // #[case::lassa_l_20("lassa_L_20")]     // TODO: missing internal node times, leaf dates not refined
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // TODO: missing internal node times, leaf dates not refined
  // #[case::rsv_a_20("rsv_a_20")]         // TODO: missing internal node times, leaf dates not refined
  // #[case::tb_20("tb_20")]               // TODO: missing internal node times, leaf dates not refined
  // #[case::zika_20("zika_20")]           // TODO: read_dates strips # from headers, name_column="#name" mismatches
  #[trace]
  #[ignore = "sparse-vs-dense discrepancy: max 0.92 years at root-adjacent node (grid-width difference)"]
  // TODO: investigate why sparse branch distribution grids produce 0.92-year shift at
  // root-adjacent nodes. All other nodes within 0.1 years. Related:
  // kb/issues/M-timetree-branch-grid-uniform-resolution.md
  fn test_gm_runner_marginal_sparse_dense_consistency(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = case.marginal_dense();

    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(case.rerooted_tree_nwk())?;

    let mut graph: GraphTimetree = graph;
    let dates = load_dates_for_dataset(dataset)?;
    let constraints = load_date_constraints(&dates, &graph, &names)?;

    let aln = load_alignment_for_dataset(dataset)?;
    let fitch = create_fitch_partition(&graph, 0, ALPHABET.clone(), &aln, &names)?;
    let sparse_partition = PartitionTimetree::Sparse(
      fitch.into_marginal_sparse(jc69(JC69Params::default())?, &graph)?,
    );

    let mut partitions: Vec<PartitionTimetree> = vec![sparse_partition];
    initialize_marginal(&graph, &profile_branch_lengths(&branch_lengths), &mut partitions, &aln, &names)?.value();
    let mut clock_state = ClockState::new(&graph);
    initialize_node_divergences(&graph, &mut clock_state, &branch_lengths, &names)?;

    let times = TimetreeState::seed_from_values(&graph, &constraints).likely_times();
    let mut clock_estimate_state = ClockState::seed_from_values(&graph, &times);
    let names_tt_1 = names.clone();
    let clock_model = estimate_clock_model_with_reroot_policy(
      &mut graph,
      &mut clock_estimate_state,
      &ClockParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      &RerootParams::default(),
      &mut branch_lengths,
      None, &names_tt_1
    )?
    .into_clock_model()?;

    let mut state = TimetreeState::new(&graph);
    let run_branch_lengths = branch_lengths;
    let run_names = names.clone();
    run_timetree(
      &mut graph,
      &mut partitions,      &run_branch_lengths,
      &run_names,
      &clock_model,
      None,
      false,
      &mut state,
      &mut clock_state,
    )?;

    let actual = extract_node_times(&graph, &names, &state);
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-6);

    Ok(())
  }
}
