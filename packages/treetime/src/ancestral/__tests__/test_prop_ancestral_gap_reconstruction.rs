#[cfg(test)]
mod tests {
  use crate::ancestral::params::MethodAncestral;
  use crate::test_utils::arb_simulated_case;
  use helpers::{GapCase, arb_gap_events};
  use proptest::prelude::{ProptestConfig, prop_assert_eq, proptest};

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn test_prop_ancestral_gap_reconstruction_parsimony_puts_each_clade_gap_below_its_branch(
      spec in arb_simulated_case(),
      events in arb_gap_events(),
    ) {
      let case = GapCase::simulate(&spec, events).expect("a gap case simulates");

      let actual = case.reconstructed_gaps(MethodAncestral::Parsimony, None).expect("parsimony runs");

      prop_assert_eq!(case.true_gaps(), actual, "gap events {:?}", case.events);
    }

    #[test]
    fn test_prop_ancestral_gap_reconstruction_dense_marginal_puts_each_clade_gap_below_its_branch(
      spec in arb_simulated_case(),
      events in arb_gap_events(),
    ) {
      let case = GapCase::simulate(&spec, events).expect("a gap case simulates");

      let actual = case.reconstructed_gaps(MethodAncestral::Marginal, Some(true)).expect("dense marginal runs");

      prop_assert_eq!(case.true_gaps(), actual, "gap events {:?}", case.events);
    }

    #[test]
    fn test_prop_ancestral_gap_reconstruction_sparse_marginal_puts_each_clade_gap_below_its_branch(
      spec in arb_simulated_case(),
      events in arb_gap_events(),
    ) {
      let case = GapCase::simulate(&spec, events).expect("a gap case simulates");

      let actual = case.reconstructed_gaps(MethodAncestral::Marginal, Some(false)).expect("sparse marginal runs");

      prop_assert_eq!(case.true_gaps(), actual, "gap events {:?}", case.events);
    }
  }

  mod helpers {
    use crate::ancestral::params::{AncestralParams, MethodAncestral};
    use crate::ancestral::pipeline;
    use crate::cancel::NoopCancel;
    use crate::gtr::get_gtr::GtrModel;
    use crate::partition::marginal::sample::SampleMode;
    use crate::progress::NoopProgress;
    use crate::seq::sink::SeqStream;
    use crate::test_utils::{CaseSpec, RecordingSeqSink, SimulatedCase, ancestral_input};
    use eyre::Report;
    use itertools::Itertools;
    use proptest::prelude::any;
    use proptest::strategy::Strategy;
    use std::collections::BTreeMap;
    use treetime_graph::node::GraphNodeKey;
    use treetime_primitives::Seq;
    use treetime_simulation::branch_lengths::exact;
    use treetime_simulation::gaps::{GAP, GapEvent, apply_gap_events, disjoint_gap_events};
    use treetime_utils::make_report;
    use treetime_utils::sync::random::get_random_number_generator;

    #[derive(Clone, Copy, Debug)]
    pub(super) struct GapEventSettings {
      seed: u64,
      count: usize,
      length: usize,
    }

    pub(super) fn arb_gap_events() -> impl Strategy<Value = GapEventSettings> {
      (any::<u64>(), 1..=6_usize, 1..=5_usize).prop_map(|(seed, count, length)| GapEventSettings {
        seed,
        count,
        length,
      })
    }

    pub(super) struct GapCase {
      case: SimulatedCase,
      sequences: BTreeMap<GraphNodeKey, Seq>,
      pub(super) events: Vec<GapEvent>,
    }

    impl GapCase {
      pub(super) fn simulate(spec: &CaseSpec, settings: GapEventSettings) -> Result<Self, Report> {
        let case = spec.simulate()?;
        let events = disjoint_gap_events(
          &case.tree.edge_durations(),
          settings.count,
          settings.length,
          spec.length,
          &mut get_random_number_generator(settings.seed),
        )?;
        let sequences = apply_gap_events(&case.tree.graph, &case.sequences.sequences, &events)?;
        Ok(Self {
          case,
          sequences,
          events,
        })
      }

      pub(super) fn true_gaps(&self) -> BTreeMap<String, Vec<(usize, usize)>> {
        self
          .case
          .tree
          .graph
          .get_internal_nodes()
          .map(|node| {
            (
              self.case.tree.names[node.key()].to_owned(),
              gap_ranges(&self.sequences[&node.key()]),
            )
          })
          .collect()
      }

      pub(super) fn reconstructed_gaps(
        &self,
        method: MethodAncestral,
        dense: Option<bool>,
      ) -> Result<BTreeMap<String, Vec<(usize, usize)>>, Report> {
        let tree = &self.case.tree;
        let params = AncestralParams {
          method,
          dense,
          include_leaves: false,
          report_ambiguous: false,
          impute_missing_data: false,
          seed: 0,
          sample_from_profile: SampleMode::Argmax,
        };
        let lengths = exact(&tree.edge_durations(), &self.case.model)?;
        let input = ancestral_input(tree, &self.sequences, &lengths, &GtrModel::Infer)?;
        let mut sink = RecordingSeqSink::default();
        pipeline::run(
          &params,
          input,
          Some(SeqStream::new(&mut sink, &tree.names)),
          &NoopCancel,
          &NoopProgress,
          &NoopProgress,
        )?;
        let emitted: BTreeMap<GraphNodeKey, Seq> = sink
          .items
          .into_iter()
          .filter(|(_, emitted, _)| *emitted)
          .map(|(node, _, seq)| (node, seq))
          .collect();
        tree
          .graph
          .get_internal_nodes()
          .map(|node| {
            let name = &tree.names[node.key()];
            let seq = emitted
              .get(&node.key())
              .ok_or_else(|| make_report!("The reconstruction emitted no sequence for node {name}"))?;
            Ok((name.to_owned(), gap_ranges(seq)))
          })
          .collect()
      }
    }

    fn gap_ranges(seq: &Seq) -> Vec<(usize, usize)> {
      seq
        .iter()
        .enumerate()
        .chunk_by(|(_, state)| **state == GAP)
        .into_iter()
        .filter_map(|(is_gap, mut sites)| {
          let (start, _) = sites.next()?;
          let end = sites.last().map_or(start, |(site, _)| site) + 1;
          is_gap.then_some((start, end))
        })
        .collect()
    }
  }
}
