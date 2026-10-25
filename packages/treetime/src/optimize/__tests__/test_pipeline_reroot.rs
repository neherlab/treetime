#[cfg(test)]
mod tests {
  use crate::cancel::NoopCancel;
  use crate::optimize::pipeline::{OptimizeInput, run};
  use crate::progress::NoopProgress;
  use crate::reroot::params::{RerootMethod, RerootSpec};
  use crate::test_utils::sequence_source;
  use eyre::Report;
  use helpers::{
    assert_branch_lengths_valid, jc69_sequences, leaf_names, load, load_with_unary_root, midpoint_root_clades,
    params_with, root_child_keys, root_clades, root_key,
  };
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::node::GraphNodeKey;

  #[test]
  fn test_optimize_pipeline_reroot_min_dev_changes_root() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let leaves_before = graph.get_leaves().count();
    let root_children_before = root_child_keys(&graph);
    let output = run(
      &params_with(Some(RerootSpec::Method(RerootMethod::MinDev))),
      OptimizeInput {
        graph,
        names,
        sequences: jc69_sequences(sequence_source(alphabet, sequences)?)?,
        branch_lengths,
      },
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().count(), leaves_before);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    let root_children_after = root_child_keys(&output.graph);
    assert_ne!(
      root_children_before, root_children_after,
      "min-dev reroot should change the root position on this unbalanced dataset"
    );
    Ok(())
  }

  #[test]
  fn test_optimize_pipeline_reroot_tips_changes_root() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let leaves_before = graph.get_leaves().count();
    let root_before = root_key(&graph);
    let tips: Vec<GraphNodeKey> = graph.get_leaves().map(|leaf| leaf.key()).take(2).collect();
    let output = run(
      &params_with(Some(RerootSpec::Tips(tips))),
      OptimizeInput {
        graph,
        names,
        sequences: jc69_sequences(sequence_source(alphabet, sequences)?)?,
        branch_lengths,
      },
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().count(), leaves_before);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    let root_after = root_key(&output.graph);
    assert_ne!(root_before, root_after, "tip-based reroot should move the root");
    Ok(())
  }

  #[test]
  fn test_optimize_pipeline_reroot_min_dev_dense_completes() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let leaves_before = graph.get_leaves().count();

    let mut params = params_with(Some(RerootSpec::Method(RerootMethod::MinDev)));
    params.dense = Some(true);
    let output = run(
      &params,
      OptimizeInput {
        graph,
        names,
        sequences: jc69_sequences(sequence_source(alphabet, sequences)?)?,
        branch_lengths,
      },
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().count(), leaves_before);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense(  true)]
  #[case::sparse( false)]
  #[trace]
  fn test_optimize_pipeline_reroot_min_dev_removes_a_unary_root(#[case] dense: bool) -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load_with_unary_root()?;
    let mut expected: Vec<String> = sequences
      .values()
      .filter(|node| node.seq.is_some())
      .map(|node| node.name.clone())
      .collect();
    expected.sort();

    let mut params = params_with(Some(RerootSpec::Method(RerootMethod::MinDev)));
    params.dense = Some(dense);
    let output = run(
      &params,
      OptimizeInput {
        graph,
        names,
        sequences: jc69_sequences(sequence_source(alphabet, sequences)?)?,
        branch_lengths,
      },
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    let mut actual = leaf_names(&output.graph, &output.names);
    actual.sort();
    assert_eq!(expected, actual);
    assert!(
      output
        .graph
        .get_nodes()
        .all(|node| &output.names[node.key()] != "UNARYROOT"),
      "the unary root is removed from the graph"
    );
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    Ok(())
  }

  #[test]
  fn test_optimize_pipeline_reroot_midpoint_roots_at_the_midpoint_of_the_starting_tree() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let expected = midpoint_root_clades()?;

    let output = run(
      &params_with(Some(RerootSpec::Method(RerootMethod::Midpoint))),
      OptimizeInput {
        graph,
        names,
        sequences: jc69_sequences(sequence_source(alphabet, sequences)?)?,
        branch_lengths,
      },
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    assert_eq!(expected, root_clades(&output.graph, &output.names)?);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    Ok(())
  }

  #[test]
  fn test_optimize_pipeline_keep_root_completes() -> Result<(), Report> {
    let (graph, names, alphabet, sequences, branch_lengths) = load()?;
    let leaves_before = graph.get_leaves().count();
    let output = run(
      &params_with(None),
      OptimizeInput {
        graph,
        names,
        sequences: jc69_sequences(sequence_source(alphabet, sequences)?)?,
        branch_lengths,
      },
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    assert_eq!(output.graph.get_leaves().count(), leaves_before);
    assert_branch_lengths_valid(&output.graph, &output.branch_lengths);
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::gtr::get_gtr::{GtrModel, JC69Params};
    use crate::optimize::indel::IndelPolicy;
    use crate::optimize::params::OptimizeLoopSettings;
    use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
    use crate::optimize::partition_rates::PartitionRates;
    use crate::optimize::pipeline::OptimizeParams;
    use crate::partition::split::SequenceInput;
    use crate::progress::NoopProgress;
    use crate::reroot::date_free::{LengthsAdvice, reroot_date_free};
    use crate::reroot::params::{RerootMethod, RerootSpec};
    use crate::seq::alignment::{NodeSeqInput, SequenceSource};
    use crate::test_utils::{TestDataset, leaf_seq_inputs};
    use eyre::Report;
    use std::collections::BTreeSet;
    use treetime_graph::clades::{TipIds, graph_clade_keys};
    use treetime_graph::key_map::{EdgeMap, NodeMap};

    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::{NwkParse, nwk_read};

    pub(super) type Loaded = (Graph, NodeNames, Alphabet, NodeMap<NodeSeqInput>, EdgeMap<Option<f64>>);

    pub(super) fn load() -> Result<Loaded, Report> {
      loaded(nwk_read(TestDataset::Flu20.tree().as_bytes())?)
    }

    pub(super) fn load_with_unary_root() -> Result<Loaded, Report> {
      let subtree = TestDataset::Flu20.tree().trim_end().trim_end_matches(';');
      loaded(nwk_read(format!("({subtree}:0.001)UNARYROOT;").as_bytes())?)
    }

    pub(super) fn jc69_sequences(source: SequenceSource) -> Result<SequenceInput, Report> {
      SequenceInput::whole_alignments(vec![source], &GtrModel::JC69(JC69Params::default()))
    }

    pub(super) fn params_with(reroot_spec: Option<RerootSpec>) -> OptimizeParams {
      OptimizeParams {
        dense: Some(false),
        initial_guess: InitialGuessMode::default(),
        reroot_spec,
        loop_settings: OptimizeLoopSettings {
          max_iter: 2,
          dp: 0.1,
          damping: 0.75,
          opt_method: BranchOptMethod::default(),
          indel_policy: IndelPolicy::Include,
          topology_ops: TopologyOps::default(),
          partition_rates: PartitionRates::Estimate,
        },
      }
    }

    pub(super) fn assert_branch_lengths_valid(graph: &Graph, branch_lengths: &EdgeMap<Option<f64>>) {
      for edge in graph.get_edges() {
        let bl = branch_lengths[&edge.key()].expect("every edge has a branch length after optimization");
        assert!(bl.is_finite() && bl >= 0.0, "invalid branch length {bl}");
      }
    }

    pub(super) fn root_key(graph: &Graph) -> GraphNodeKey {
      graph.get_exactly_one_root().unwrap().key()
    }

    pub(super) fn root_child_keys(graph: &Graph) -> Vec<GraphNodeKey> {
      let root = graph.get_exactly_one_root().unwrap();
      root
        .outbound()
        .iter()
        .map(|&edge_key| graph.get_edge(edge_key).unwrap().target())
        .collect()
    }

    pub(super) fn leaf_names(graph: &Graph, names: &NodeNames) -> Vec<String> {
      graph.get_leaves().map(|leaf| names[leaf.key()].to_owned()).collect()
    }

    pub(super) fn midpoint_root_clades() -> Result<BTreeSet<BTreeSet<String>>, Report> {
      let (mut graph, mut names, _, _, mut branch_lengths) = load()?;
      reroot_date_free(
        &mut graph,
        &mut branch_lengths,
        &RerootSpec::Method(RerootMethod::Midpoint),
        &|_| false,
        LengthsAdvice::InitialBranchLengths,
        &mut names,
        &NoopProgress,
      )?;
      root_clades(&graph, &names)
    }

    pub(super) fn root_clades(graph: &Graph, names: &NodeNames) -> Result<BTreeSet<BTreeSet<String>>, Report> {
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

    fn loaded(nwk_parsed: NwkParse) -> Result<Loaded, Report> {
      let alphabet = Alphabet::default();
      let names = nwk_parsed.names.clone();
      let sequences = leaf_seq_inputs(&nwk_parsed.graph, &names, TestDataset::Flu20.alignment()?);
      Ok((nwk_parsed.graph, names, alphabet, sequences, nwk_parsed.branch_lengths))
    }
  }
}
