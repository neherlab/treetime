#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::cancel::NoopCancel;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::{ClockFit, ClockVarianceParams};
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::partition::create::{Representation, build_marginal_partition};
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::test_utils::{marginal_timetree_params, point_date_constraints};
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::params::{TimeMarginalMode, TimetreeContext};
  use crate::timetree::pre_loop::{PreLoopInputs, PreLoopState, run_pre_loop};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  const INPUT_BRANCH_LENGTH: f64 = 0.5;

  const SEQUENCE_LENGTH: f64 = 24.0;

  #[test]
  fn test_pre_loop_ml_step_shortens_branches_far_longer_than_the_alignment_supports() -> Result<(), Report> {
    let length = INPUT_BRANCH_LENGTH;
    let newick = format!("((A:{length},B:{length})AB:{length},(C:{length},D:{length})CD:{length})root;");
    let nwk_parsed = nwk_read_str(&newick)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
      indoc! {"
        >A
        ACGTACGTACGTACGTACGTACGT
        >B
        ACGTACGTACGAACGTACGTACGT
        >C
        ACGTTCGTACGTACGTACGTACGT
        >D
        ACGTTCGTACGTACGTACGTACGT
      "},
      &alphabet,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
    let reconstruction = build_marginal_partition(
      Representation::Sparse,
      GtrModelName::JC69,
      &graph,
      0,
      alphabet,
      &node_seq_inputs(&graph, &names, aln),
      &branch_lengths_or_zero(&branch_lengths),
      &NoopProgress,
    )?;
    let date_constraints = point_date_constraints(&graph, &names, &[("A", 2010.0), ("B", 2011.0), ("C", 2012.0)]);
    let context = TimetreeContext {
      time_marginal: TimeMarginalMode::Never,
      date_constraints,
      covariation_clock_params: ClockVarianceParams::default(),
      branch_params: BranchPointOptimizationParams::default(),
    };
    let params = marginal_timetree_params();
    let inputs = PreLoopInputs {
      params: &params,
      context: &context,
      names: &names,
      has_alignment: true,
    };
    let clock_fit = ClockFit {
      model: ClockModel::for_testing(0.001, 0.0),
      points: vec![],
    };
    let input_edges: BTreeSet<GraphEdgeKey> = branch_lengths.keys().copied().collect();
    let state = PreLoopState::new(graph, branch_lengths, BranchModel::Marginal(reconstruction), clock_fit);

    let state = run_pre_loop(&inputs, state, &NoopCancel, &NoopProgress, &NoopProgress)?;

    assert_eq!(
      input_edges,
      state.branch_lengths.keys().copied().collect::<BTreeSet<_>>()
    );
    let not_shortened: BTreeMap<GraphEdgeKey, Option<f64>> = state
      .branch_lengths
      .iter()
      .filter(|(_, length)| !length.is_some_and(|length| length < INPUT_BRANCH_LENGTH))
      .map(|(key, length)| (*key, *length))
      .collect();
    assert_eq!(
      BTreeMap::new(),
      not_shortened,
      "with at most 2 differences in {SEQUENCE_LENGTH} sites, no branch can stay at {INPUT_BRANCH_LENGTH}"
    );
    assert_eq!(BTreeSet::new(), state.outliers);
    assert_eq!(None, state.filter_divergences);
    Ok(())
  }
}
