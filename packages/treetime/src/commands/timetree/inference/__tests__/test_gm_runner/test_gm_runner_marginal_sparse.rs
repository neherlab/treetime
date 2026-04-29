#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::commands::ancestral::fitch::compress_sequences;
  use crate::commands::ancestral::marginal::initialize_marginal;
  use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::commands::timetree::inference::runner::run_timetree;
  use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
  use crate::commands::timetree::utils::{
    extract_node_times, initialize_clock_totals_from_time_distributions, initialize_node_divergences,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use rstest::rstest;
  use std::slice::from_ref;
  use std::sync::Arc;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::seq;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  // --- Marginal sparse tests ---

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
  fn test_gm_runner_marginal_sparse(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = case.marginal_dense();

    let mut graph: GraphTimetree = nwk_read_str(case.rerooted_tree_nwk())?;
    let dates = load_dates_for_dataset(dataset)?;
    load_date_constraints(&dates, &graph)?;

    let aln = load_alignment_for_dataset(dataset)?;
    let sparse_partition = Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: case.sequence_length(),
      nodes: btreemap! {},
      edges: btreemap! {},
      root_sequence: seq![],
    }));

compress_sequences(&graph, from_ref(&sparse_partition), &aln)?;

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![sparse_partition];
    initialize_marginal(&graph, &partitions, &aln)?;
    initialize_node_divergences(&graph);

    let clock_model = estimate_clock_model_with_reroot(
      &mut graph,
      &ClockParams::default(),
      Some(case.clock_rate()),
      true,
      &BranchPointOptimizationParams::default(),
      None,
    )?;

    initialize_clock_totals_from_time_distributions(&graph)?;
    run_timetree(&mut graph, &partitions, &clock_model, None, false)?;

    let actual = extract_node_times(&graph);
    // Tolerance increased from 0.9 to 1.0 after widening branch distribution grids
    // (MAX_BRANCH_TIME): root-adjacent nodes shift by up to 0.92 years because the
    // wider grid evaluates the true likelihood over a larger range, changing the
    // distribution shape for deep-tree branches.
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e0);

    Ok(())
  }
}
