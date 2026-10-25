#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::cancel::NoopCancel;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_regression::ClockFit;
  use crate::clock::gaussian::GaussianTree;
  use crate::gtr::get_gtr::{GtrModel, JC69Params};
  use crate::partition::create::{Representation, build_marginal_partition};
  use crate::partition::inferred_rate::InferredRate;
  use crate::partition::set::PartitionSet;
  use crate::pretty_assert_ulps_eq;
  use crate::progress::NoopProgress;
  use crate::test_utils::leaf_seq_inputs;
  use crate::test_utils::{jc69_input_branch_model, marginal_timetree_params, point_date_constraints};
  use crate::timetree::branch_model::BranchModel;
  use crate::timetree::inference::bad_branches::bad_leaves;
  use crate::timetree::inference::result::given_times;
  use crate::timetree::params::{FinalSequences, TimetreeContext, TimetreeParams};
  use crate::timetree::pre_loop::{PreLoopInputs, PreLoopState, run_pre_loop};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  const INPUT_BRANCH_LENGTH: f64 = 0.5;

  const SEQUENCE_LENGTH: f64 = 24.0;

  #[test]
  fn test_pre_loop_ml_step_shortens_branches_far_longer_than_the_alignment_supports() -> Result<(), Report> {
    let length = INPUT_BRANCH_LENGTH;
    let newick = format!("((A:{length},B:{length})AB:{length},(C:{length},D:{length})CD:{length})root;");
    let nwk_parsed = nwk_read(newick.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln: Vec<AlignmentRecord> = fasta_read(
      indoc! {"
        >A
        ACGTACGTACGTACGTACGTACGT
        >B
        ACGTACGTACGAACGTACGTACGT
        >C
        ACGTTCGTACGTACGTACGTACGT
        >D
        ACGTTCGTACGTACGTACGTACGT
      "}
      .as_bytes(),
      &alphabet,
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();
    let reconstruction = PartitionSet::from(build_marginal_partition(
      Representation::Sparse,
      &GtrModel::JC69(JC69Params::default()),
      InferredRate::Unit,
      &graph,
      alphabet,
      leaf_seq_inputs(&graph, &names, aln),
      &branch_lengths_or_zero(&branch_lengths),
      &NoopProgress,
    )?);
    let date_constraints = point_date_constraints(&graph, &names, &[("A", 2010.0), ("B", 2011.0), ("C", 2012.0)]);
    let context = TimetreeContext {
      final_sequences: FinalSequences::Skip,
      date_constraints,
      sequence_length: 16,
    };
    let params = marginal_timetree_params();
    let inputs = PreLoopInputs {
      params: &params,
      context: &context,
      has_alignment: true,
    };
    let clock_fit = ClockFit {
      model: ClockModel::for_testing(0.001, 0.0),
      points: vec![],
    };
    let input_edges: BTreeSet<GraphEdgeKey> = branch_lengths.keys().copied().collect();
    let state = PreLoopState::new(
      graph,
      names,
      branch_lengths,
      BranchModel::Marginal(reconstruction),
      clock_fit,
    );

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
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::outgroup_slope_below_the_floor(helpers::Root::Outgroup, helpers::SLOPE_BELOW_THE_FLOOR)]
  #[case::oldest_slope_below_the_floor(  helpers::Root::Oldest,   helpers::SLOPE_BELOW_THE_FLOOR)]
  #[case::oldest_negative_slope(         helpers::Root::Oldest,   helpers::NEGATIVE_SLOPE)]
  #[trace]
  fn test_pre_loop_refit_on_a_root_without_slope_search_falls_back_to_the_gaussian_rate(
    #[case] root: helpers::Root,
    #[case] newick: &str,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(newick.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let date_constraints = point_date_constraints(&graph, &names, &helpers::TIP_DATES);
    let context = TimetreeContext {
      final_sequences: FinalSequences::Skip,
      date_constraints,
      sequence_length: helpers::SEQUENCE_LENGTH,
    };
    let params = TimetreeParams {
      reroot_spec: Some(root.spec(&graph, &names)),
      ..marginal_timetree_params()
    };
    let inputs = PreLoopInputs {
      params: &params,
      context: &context,
      has_alignment: false,
    };
    let clock_fit = ClockFit {
      model: ClockModel::for_testing(0.001, 0.0),
      points: vec![],
    };
    let state = PreLoopState::new(
      graph,
      names,
      nwk_parsed.branch_lengths,
      jc69_input_branch_model(helpers::SEQUENCE_LENGTH),
      clock_fit,
    );

    let state = run_pre_loop(&inputs, state, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let tip_dates = given_times(&state.graph, &context.date_constraints)?;
    let excluded = bad_leaves(&state.graph, &context.date_constraints, &BTreeSet::new());
    let gaussian_rate = GaussianTree::new(
      &state.graph,
      &state.branch_lengths,
      helpers::SEQUENCE_LENGTH,
      &tip_dates,
      &excluded,
    )?
    .rate()?;
    pretty_assert_ulps_eq!(gaussian_rate, state.clock_fit.model.clock_rate(), max_ulps = 0);
    Ok(())
  }

  mod helpers {
    use crate::reroot::params::{RerootMethod, RerootSpec};
    use crate::test_utils::find_node_key_by_name;
    use treetime_graph::graph::Graph;
    use treetime_graph::node_names::NodeNames;

    pub(super) const SEQUENCE_LENGTH: usize = 1000;

    pub(super) const TIP_DATES: [(&str, f64); 6] = [
      ("X", 2000.0),
      ("A", 2002.0),
      ("B", 2004.0),
      ("C", 2006.0),
      ("D", 2008.0),
      ("E", 2010.0),
    ];

    pub(super) const SLOPE_BELOW_THE_FLOOR: &str = "(X:0.1,(A:0.001,B:0.001,C:0.0012,D:0.001,E:0.0014)S:0.1)root;";

    pub(super) const NEGATIVE_SLOPE: &str = "(X:0.1,(A:0.01,B:0.008,C:0.006,D:0.004,E:0.002)S:0.1)root;";

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Root {
      Outgroup,
      Oldest,
    }

    impl Root {
      pub(super) fn spec(self, graph: &Graph, names: &NodeNames) -> RerootSpec {
        match self {
          Self::Outgroup => RerootSpec::Tips(vec![find_node_key_by_name(graph, names, "X").expect("outgroup tip")]),
          Self::Oldest => RerootSpec::Method(RerootMethod::Oldest),
        }
      }
    }
  }
}
