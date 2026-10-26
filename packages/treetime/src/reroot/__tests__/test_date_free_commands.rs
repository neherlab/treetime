#[cfg(test)]
mod tests {
  use eyre::Report;
  use helpers::{Rooted, Spec, clock_root, date_free_root, optimize_root, timetree_root};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::min_dev( Spec::MinDev)]
  #[case::outgroup(Spec::Outgroup)]
  #[trace]
  fn test_date_free_commands_root_once_at_the_root_of_the_search_on_the_starting_tree(
    #[case] spec: Spec,
  ) -> Result<(), Report> {
    let expected = Rooted {
      clades: date_free_root(spec)?,
      reroots: 1,
    };

    let actual = (clock_root(spec)?, timetree_root(spec)?, optimize_root(spec)?);

    assert_eq!((expected.clone(), expected.clone(), expected), actual);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::cancel::NoopCancel;
    use crate::clock::pipeline::{ClockInput, ClockParams};
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, OptimizeLoopSettings, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::pipeline::{OptimizeInput, OptimizeParams};
    use crate::partition::split::SequenceInput;
    use crate::progress::{LogLevel, NoopProgress};
    use crate::reroot::date_free::{LengthsAdvice, reroot_date_free};
    use crate::reroot::params::{RerootMethod, RerootSpec};
    use crate::test_utils::{
      RecordingLog, TestDataset, dates_by_node, find_node_key_by_name, leaf_seq_inputs, marginal_timetree_params,
      sequence_source,
    };
    use crate::timetree::params::TimetreeParams;
    use crate::timetree::pipeline::{TimetreeInput, TimetreeSequences};
    use crate::{clock, optimize, timetree};
    use eyre::Report;
    use std::collections::BTreeSet;
    use treetime_graph::clades::{TipIds, graph_clade_keys};
    use treetime_graph::graph::Graph;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkParse, nwk_read};

    const DATASET: TestDataset = TestDataset::Zika20;

    const CLOCK_FILTER_IQD: f64 = 3.0;

    const OUTGROUP: [&str; 2] = [
      "BeH818995|KU365777|2015-07-01|brazil",
      "SSABR1|KU707826|2015-07-01|brazil",
    ];

    pub(super) type RootClades = BTreeSet<BTreeSet<String>>;

    #[derive(Clone, Debug, PartialEq, Eq)]
    pub(super) struct Rooted {
      pub(super) clades: RootClades,
      pub(super) reroots: usize,
    }

    #[derive(Clone, Copy, Debug)]
    pub(super) enum Spec {
      MinDev,
      Outgroup,
    }

    impl Spec {
      fn reroot_spec(self, tree: &NwkParse) -> RerootSpec {
        match self {
          Self::MinDev => RerootSpec::Method(RerootMethod::MinDev),
          Self::Outgroup => RerootSpec::Tips(
            OUTGROUP
              .iter()
              .map(|tip| find_node_key_by_name(&tree.graph, &tree.names, tip).expect("the outgroup tip exists"))
              .collect(),
          ),
        }
      }
    }

    pub(super) fn date_free_root(spec: Spec) -> Result<RootClades, Report> {
      let tree = nwk_read(DATASET.tree().as_bytes())?;
      let reroot_spec = spec.reroot_spec(&tree);
      let NwkParse {
        mut graph,
        mut names,
        mut branch_lengths,
        ..
      } = tree;
      reroot_date_free(
        &mut graph,
        &mut branch_lengths,
        &reroot_spec,
        &|_| false,
        LengthsAdvice::OptimizeTheTree,
        &mut names,
        &NoopProgress,
      )?;
      root_clades(&graph, &names)
    }

    pub(super) fn clock_root(spec: Spec) -> Result<Rooted, Report> {
      let log = RecordingLog::default();
      let tree = nwk_read(DATASET.tree().as_bytes())?;
      let params = ClockParams {
        clock_filter: CLOCK_FILTER_IQD,
        allow_negative_rate: false,
        reroot_spec: Some(spec.reroot_spec(&tree)),
        prune_filtered_tips: None,
      };
      let input = ClockInput {
        dates: dates_by_node(DATASET.dates()?, &tree.graph, &tree.names),
        graph: tree.graph,
        names: tree.names,
        branch_lengths: tree.branch_lengths,
      };
      let output = clock::pipeline::run(&params, input, &NoopCancel, &NoopProgress, &log)?;
      rooted(&output.graph, &output.names, &log)
    }

    pub(super) fn timetree_root(spec: Spec) -> Result<Rooted, Report> {
      let log = RecordingLog::default();
      let tree = nwk_read(DATASET.tree().as_bytes())?;
      let params = TimetreeParams {
        reroot_spec: Some(spec.reroot_spec(&tree)),
        clock_filter: CLOCK_FILTER_IQD,
        prune_filtered_tips: None,
        ..marginal_timetree_params()
      };
      let input = TimetreeInput {
        dates: Some(dates_by_node(DATASET.dates()?, &tree.graph, &tree.names)),
        sequences: TimetreeSequences::Alignments(jc69_alignment(&tree)?),
        graph: tree.graph,
        names: tree.names,
        branch_lengths: tree.branch_lengths,
      };
      let output = timetree::pipeline::run(&params, input, None, None, &NoopCancel, &NoopProgress, &log)?;
      rooted(&output.graph, &output.names, &log)
    }

    pub(super) fn optimize_root(spec: Spec) -> Result<Rooted, Report> {
      let log = RecordingLog::default();
      let tree = nwk_read(DATASET.tree().as_bytes())?;
      let params = OptimizeParams {
        dense: Some(false),
        initial_guess: InitialGuessMode::default(),
        reroot_spec: Some(spec.reroot_spec(&tree)),
        loop_settings: OptimizeLoopSettings {
          max_iter: 2,
          dp: 0.1,
          damping: 0.75,
          opt_method: BranchOptMethod::default(),
          indel_policy: IndelPolicy::Include,
          topology_ops: TopologyOps::default(),
          partition_rates: PartitionRates::Estimate,
        },
      };
      let input = OptimizeInput {
        sequences: jc69_alignment(&tree)?,
        graph: tree.graph,
        names: tree.names,
        branch_lengths: tree.branch_lengths,
      };
      let output = optimize::pipeline::run(&params, input, &NoopCancel, &NoopProgress, &log)?;
      rooted(&output.graph, &output.names, &log)
    }

    fn jc69_alignment(tree: &NwkParse) -> Result<SequenceInput, Report> {
      let records = leaf_seq_inputs(&tree.graph, &tree.names, DATASET.alignment()?);
      let source = sequence_source(Alphabet::default(), records)?;
      SequenceInput::whole_alignments(vec![source], &GtrModel::JC69(JC69Params::default()))
    }

    fn rooted(graph: &Graph, names: &NodeNames, log: &RecordingLog) -> Result<Rooted, Report> {
      let reroots = log
        .messages(LogLevel::Info)
        .iter()
        .filter(|message| message.starts_with("Rerooted to "))
        .count();
      Ok(Rooted {
        clades: root_clades(graph, names)?,
        reroots,
      })
    }

    fn root_clades(graph: &Graph, names: &NodeNames) -> Result<RootClades, Report> {
      let mut ids = TipIds::default();
      let clades = graph_clade_keys(graph, names, &mut ids)?;
      let schedule = GraphSchedule::new(graph)?;
      let root_key = graph.get_exactly_one_root()?.key();
      Ok(
        schedule
          .children(root_key)
          .iter()
          .map(|(child_key, _)| clades[child_key].iter().map(|&id| ids.name(id).to_owned()).collect())
          .collect(),
      )
    }
  }
}
