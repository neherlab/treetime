#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{
    ALPHABET, OUTPUTS, load_alignment_for_dataset, load_dates_for_dataset,
  };
  use crate::commands::ancestral::marginal::initialize_marginal;
  use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot};
  use crate::commands::clock::date_constraints::load_date_constraints;
  use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::commands::timetree::inference::runner::run_timetree;
  use crate::representation::partition::traits::PartitionTimetreeAll;
  use crate::commands::timetree::utils::{
    extract_node_times, initialize_clock_totals_from_time_distributions, initialize_node_divergences,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::timetree::GraphTimetree;
  use crate::representation::payload::timetree::{EdgeTimetree, NodeTimetree};
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  // --- Marginal dense tests ---

  #[rustfmt::skip]
#[rstest]
  // #[case::dengue_20("dengue_20")]       // TODO: missing internal node times, leaf dates not refined
  // #[case::ebola_20("ebola_20")]         // TODO: golden master node key mismatch (v0 captures 11 internal nodes, v1 rerooting produces 19)
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // #[case::lassa_l_20("lassa_L_20")]     // TODO: missing internal node times, leaf dates not refined
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // TODO: missing internal node times, leaf dates not refined
  // #[case::rsv_a_20("rsv_a_20")]         // TODO: missing internal node times, leaf dates not refined
  // #[case::tb_20("tb_20")]               // TODO: missing internal node times, leaf dates not refined
  // #[case::zika_20("zika_20")]           // TODO: read_dates strips # from headers, name_column="#name" mismatches
  #[trace]
  #[ignore = "dense-vs-v0 discrepancy: max 0.92 years at root-adjacent node (grid-width difference)"]
  // TODO: investigate why uniform branch distribution grids produce 0.92-year shift at
  // root-adjacent nodes. All other nodes within 0.1 years. Related:
  // kb/issues/M-timetree-branch-grid-uniform-resolution.md
  fn test_gm_runner_marginal_dense(#[case] dataset: &str) -> Result<(), Report> {
    let case = &OUTPUTS[dataset];
    let expected = case.marginal_dense();

    let mut graph: GraphTimetree = nwk_read_str(case.rerooted_tree_nwk())?;
    let dates = load_dates_for_dataset(dataset)?;
    load_date_constraints(&dates, &graph)?;

    let aln = load_alignment_for_dataset(dataset)?;
    let dense_partition = Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: ALPHABET.clone(),
      length: case.sequence_length(),
      nodes: btreemap! {},
      edges: btreemap! {},
    }));

    let partitions: Vec<Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>> = vec![dense_partition];
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
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-6);

    Ok(())
  }
}
